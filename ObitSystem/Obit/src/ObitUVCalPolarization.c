/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2016                                          */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitUVDesc.h"
#include "ObitUVCalPolarizationDef.h"
#include "ObitUVCalPolarization.h"
#include "ObitTablePD.h"
#include "ObitSinCos.h"
#include "ObitComplex.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVCalPolarization.c
 * ObitUVCal utilities for applying polarization calibration to uv data.
 */

/*-------------------Private function prototypes-------------------*/
/** Private:  Create structure for polarization calibration . */
static ObitUVCalPolarizationS* newObitUVCalPolarizationS (ObitUVCal *in);

/** Private: Update calibration arrays. */
static void 
ObitUVCalPolarizationUpdate (ObitUVCalPolarizationS *in, ObitUVCal *UVCal, 
			     ObitUVCalCalibrateS *cal, 
			     ofloat time, olong SourID, olong SubA, olong FreqID,
			     ObitErr *err);

/** Private: Set baseline IF inverse Mueller matrix for R/L Linear model */
static void LinPolIF(ObitUVCalPolarizationS *in, ObitUVCal *UVCal, olong SubA,
		     olong iant1, olong iant2, olong iChan, ObitErr *err);

/** Private: Set baseline Channel/IF inverse Mueller matrix for R/L Linear model */
static void LinPolCh(ObitUVCalPolarizationS *in, ObitUVCal *UVCal, 
		     olong iant1, olong iant2, olong iChan, ObitErr *err);

/** Private: Set baseline inverse Mueller matrix for Elipticity/Orientation model */
static void OriPol(ObitUVCalPolarizationS *in, olong iant1, olong iant2, olong iChan, ObitErr *err);

/** Private: Set baseline inverse Mueller matrix for R/L Linear model resolved source */
static void VLBIPol(ObitUVCalPolarizationS *in, olong iant1, olong iant2, olong iChan, ObitErr *err);

/** Private: Set baseline inverse Mueller matrix for X/Y Linear model  */
static void LXYPol(ObitUVCalPolarizationS *in, olong iant1, olong iant2, olong iChan, ObitErr *err);

/** Private: Set inverse antenna Jones matrix per IF */
static void SetInvJonesIF(ObitUVCalPolarizationS *in, ObitAntennaList *Ant,
			  ObitUVCalCalibrateS *cal, olong iChan, olong iant);

/** Private: Set inverse antenna Jones matrix per channel/IF */
static void SetInvJonesCh(ObitUVCalPolarizationS *in, ObitUVCalCalibrateS *cal, 
			  olong iChan, olong iant, ObitErr *err);
/** Private: 4x4 Matrix * 4x1 vector multiply */
static void MatxVec4Mult(ofloat* in1, ofloat* in2, ofloat* out);
/** Private: Muller matrix from outer product of Jones matrices */
static void MatxOuter(ofloat* in1, ofloat* in2, ofloat* out);
/*----------------------Public functions---------------------------*/
/**
 * Initialize structures for polarization Calibration .
 * \param in   Calibration Object.
 * \param sel  Data selector.
 * \param desc Data descriptor.
 * \param err  ObitError stack.
 */
void ObitUVCalPolarizationInit (ObitUVCal *in, ObitUVSel *sel, ObitUVDesc *desc, 
				ObitErr *err)
{
  /*  ObitIOCode retCode;*/
  ObitUVCalPolarizationS *me;
  ofloat t1=0.0, t2, t3;
  olong i, size;
  gchar *routine = "ObitUVCalPolarizationInit";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVCalIsA(in));

  in->doPol = sel->doPolCal;
  if (!in->doPol) return;

  in->polnCal = newObitUVCalPolarizationS(in);

  /* pointer to calibration structure */
  me = in->polnCal;

  /* Copy Selector information */
  me->bChan       = sel->startChann;
  me->eChan       = sel->startChann + sel->numberChann - 1;
  me->bIF         = sel->startIF;
  me->eIF         = sel->startIF + sel->numberIF - 1;
  me->SubA        = sel->SubA;
  me->FreqID      = sel->FreqID;

  /* Copy Cal information */
  me->numSubA     = in->numSubA;

  /* Copy descriptor information */
  me->numAnt    = desc->maxAnt;
  me->numSubA   = desc->numSubA;
  if (desc->jlocif>=0) me->numIF = desc->inaxes[desc->jlocif];
  else                 me->numIF = 1;
  me->numChan   = desc->inaxes[desc->jlocf];
  me->polType   = OBIT_UVPoln_Unknown;

  /* Are all stokes requires - no for now */
  me->allStokes = FALSE;

  /* Create structures */
  /* Poln cal structure from AIPS PD table */
  me->perChan = in->PDVer>=0;
  if (me->perChan) {
    me->PCal = ObitPolCalListCreate("PolCal", in->PDTable, err);
    if (err->error) Obit_traceback_msg (err, routine, me->name);
    me->polType = me->PCal->polType;  /* Polarization parameterization type */
  }

  me->curPA    = g_realloc(me->curPA,    me->numAnt*sizeof(ofloat));
  me->curCosPA = g_realloc(me->curCosPA, me->numAnt*sizeof(ofloat));
  me->curSinPA = g_realloc(me->curSinPA, me->numAnt*sizeof(ofloat));
  me->Jones    = g_realloc(me->Jones,    me->numAnt*sizeof(ofloat*));
  /* Entry per channel or IF? */
  if (me->perChan) size = 8 * me->numIF * me->numChan;
  else             size = 8 * me->numIF;
  for (i=0; i<me->numAnt; i++) me->Jones[i] = g_malloc0(size*sizeof(ofloat));
  /* Entry per channel or IF? */
  if (me->perChan) size = 32 * me->numIF * me->numChan;
  else             size = 32 * me->numIF;
  me->PolCal   = g_realloc(me->PolCal, size*sizeof(ofloat));

  /* Init time, source */
  me->curTime   = -1.0e20;
  me->curSourID = -1;
  me->curSubA   = -1;

  /* Init Sine/cosine function */
  ObitSinCosCalc(t1, &t2, &t3);

  /* Are the data from circular or linear feeds */
  me->circFeed = desc->crval[desc->jlocs] > -4.0;

} /*  end ObitUVCalPolarizationInit */

/**
 * Polarization calibrate data.
 * For OBIT_UVPoln_VLBI and OBIT_UVPoln_ELORI, the parallactic angle correction 
 * is assumed to have already been applied.
 * Polarization calibration structure, values in order:
 *   By IF (EIF-BIF+1)
 *      A 4x4 complex matrix to be multiplied by  the observed polarization vector
 *      (RR,LL,RL,LR or XX,YY,XY,YX) to produce the corrected data.
 * Adapted from the AIPSish DATPOL.FOR
 *
 * Notes on parallactic angle:
 * For VLA style calibration (OBIT_UVPoln_Approx), no correction is made to the data 
 * early on to remove the parallactic angle.  In this case, the instrumental polarization 
 * is constant, i.e.correction terms are not a function of parallactic angle whereas the 
 * corrected data must have RL,LR rotated by the parallactic angle.
 * The VLBI stype calibration (OBIT_UVPoln_ELORI, OBIT_UVPoln_VLBI) require that the data 
 * be corrected for the effects of parallactic angle prior to any self calibration.
 * In this case, the polarization corrections rotate with parallactic angle but the 
 * corrected data need no further rotation.
 * \param in    Polarization Calibration Object.
 * \param time  Time of datum
 * \param ant1  first antenna number of baseline
 * \param ant2  second antanna of baseline.
 * \param RP    Random parameters array.
 * \param visIn 1 visibility as an array of floats
 * \param err   ObitError stack.
 */
void ObitUVCalPolarization (ObitUVCal *in, float time, olong ant1, olong ant2, 
			    ofloat *RP, ofloat *visIn, ObitErr *err)
{
  olong SubA, SourID, FreqID, iChan, limit, index, jndex, loff;
  olong i, iif, ifreq, nch, ipol, ioff, joff, koff, jrl, jlr, ia1, ia2, ifoff;
  olong j, voff[4], choff, chdelta, it1, it2;
  gboolean wflag, someOK;
  ofloat ytemp[8], dtemp[8], Lambda2, fblank = ObitMagicF();
  ofloat gr, gi, gr1, gi1, tr, ti;
  ObitUVCalPolarizationS *me;
  ObitUVCalCalibrateS   *cal;
  ObitUVDesc *desc;
  ObitUVSel *sel;
  gchar *routine="ObitUVCalPolarization";

  if (err->error) return;

  /* local pointers for structures */
  me   = in->polnCal;
  cal  = in->ampPhaseCal;
  desc = in->myDesc;
  sel  = in->mySel;

  /* Make sure some good data */
  someOK = FALSE;
  for (i=0; i<desc->ncorr; i++) someOK = someOK || (visIn[i*3+2]>0.0);
  if (!someOK) return;

  /* Subarray number in data */
  ObitUVDescGetAnts(desc, RP, &it1, &it2, &SubA);
  SubA = MIN (SubA, in->numANTable);
  nch  = in->eChan - in->bChan + 1;
  ia1  = ant1 - 1;
  ia2  = ant2 - 1;

  /* Make sure have polarization parameterization */
  if ((me->polType==OBIT_UVPoln_Unknown)  || (me->polType==OBIT_UVPoln_NoCal))
    me->polType = in->antennaLists[SubA-1]->polType;

   /* Data Freq id */
  if (desc->ilocfq >= 0) FreqID = RP[desc->ilocfq] + 0.1;
  else  FreqID = 0;

  /* Source ID */
  if (desc->ilocsu >= 0) SourID = RP[desc->ilocsu] + 0.1;
  else SourID = 0;

  /* Time for new parallactic angles? Update every 10 seconds. */
  if ((time > (me->curTime+10.0/86400.0)) || (SourID!=me->curSourID) || (SubA!=me->curSubA))
    ObitUVCalPolarizationUpdate(me, in, cal, time, SourID, SubA, FreqID, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Use central frequency channel, each IF */
  iChan = 0.5 * (in->bChan + in->eChan);

  /* Get inverse Mueller matrix by type. */
  switch (me->polType) {
  case OBIT_UVPoln_ELORI:   /* Elipticity-orientation  */
    OriPol (me, ant1, ant2, iChan, err);
    break;
  case OBIT_UVPoln_Approx:  /* R/L Linear D-term approximation */
    if (me->perChan) LinPolCh (me, in, ant1, ant2, iChan, err);
    else             LinPolIF (me, in, SubA, ant1, ant2, iChan, err);
    break;
  case OBIT_UVPoln_VLBI:    /* R/L Linear D-term approximation for resolved sources */
    VLBIPol (me, ant1, ant2, iChan, err);
    break;
  case OBIT_UVPoln_XYLin:   /* X/Y Linear D-term approximation  */
    LXYPol (me, ant1, ant2, iChan, err);
    break;
  default:
    g_assert_not_reached(); /* unknown, barf */
  }; /* end switch */
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Set order of data for matrix multiply - the matrices are for RR,RL,LR,LL
     and the data is RR,LL,RL,LR */
  voff[0] = 0;
  voff[1] = 2*desc->incs;
  voff[2] = 3*desc->incs;
  voff[3] = desc->incs;

  /* Calibration per IF or channel? */
  if (me->perChan) chdelta = me->eChan - me->bChan + 1; /* per channel/IF */
  else             chdelta = 1;                         /* per IF */

  /* loop thru if */
  ifoff = 0; /* Offset to beginning of IF matrix in PolCal */
  for (iif=  me->bIF; iif<=me->eIF; iif++) { /* loop 400 */
    ioff = (iif-1) * desc->incif;

    /* loop thru channels */
    choff = 0; /* Offset to beginning of channel matrix in PolCal for perChan */
    for (ifreq=me->bChan; ifreq<=me->eChan; ifreq++) {  /* loop 300 */
      joff = ((ifreq-1) * desc->incf + ioff); /* Offset of RR (or XX) */


      /* deal with case of missing  parallel poln; use one present for correction.
	 The "faked" data will still be flagged in the output data and is merely used
	 to correct the valid data */
      wflag = FALSE;

      koff = joff + desc->incs; /* offset to LL (or YY)*/
      if (visIn[joff+2] <= 0.0) {  /* 1st par. poln missing - use second */	
	visIn[joff]   = visIn[koff];
	visIn[joff+1] = visIn[koff+1];
	visIn[joff+2] = 0.0;
      } 

      if (visIn[koff+2] <= 0.0) { /* 2nd par. poln missing */
	visIn[koff]   = visIn[joff];
	visIn[koff+1] = visIn[joff+1];
	visIn[koff+2] = 0.0;

	/* flag all if neither parallel poln present. */
	if (visIn[joff+2] <= 0.0) wflag = TRUE;
      }

      /* Check for missing cross-hand data */
      if (in->numStok > 2) {
	jrl = joff + 2*desc->incs;
	jlr = joff + 3*desc->incs;
	
	if ((visIn[jrl+2] <= 0.0) && (visIn[jlr+2] <= 0.0)) { /* both missing */
	  /* zero cross-hand data used in correction */
	  visIn[jrl] = 0.0;  visIn[jrl+1] = 0.0;
	  visIn[jlr] = 0.0;  visIn[jlr+1] = 0.0;
	  
	} else if (visIn[jrl+2] <= 0.0) { /* 1st cross-hand missing */
	  /* Use rl=conjg(lr) approx. */
	  visIn[jrl]   =  visIn[jlr];
	  visIn[jrl+1] = -visIn[jlr+1];
	  
	} else if (visIn[jlr+2] <= 0.0) { /* 2nd cross-hand missing */
	  /* use lr=conjg(rl) approx. */
	  visIn[jlr]   =  visIn[jrl];
	  visIn[jlr+1] = -visIn[jrl+1];
	} 
      } /* end check for missing crosshand data */
      
      /* if allStokes check for any missing correlations */
      if (me->allStokes) {
	limit = in->numStok;
	jndex = joff + 2;
	for (ipol= 0; ipol< limit; ipol++) {
	  if (visIn[jndex] <= 0.0) wflag = TRUE;
	  jndex += desc->incs;
	}
      } 

      /* check for blanked IFR */
      if ((cal!=NULL) && ((cal->IFR[ia1] == fblank)  ||  (cal->IFR[ia2] == fblank)))
	wflag = TRUE;
     
      /* flag all output data if both par. hands missing or (dopol > 2) and any polzn. 
	 correlations are missing or if IFR corrections are blanked */
      
      if (wflag) {
	limit = in->numStok;
	jndex = joff + 2;
	for (ipol= 0; ipol< limit; ipol++) { /* loop 120 */
	  visIn[jndex] = 0.0;
	  jndex += desc->incs;
	} /* end loop L120:  */;
      } 
      
      /* Save in reordered array */
      dtemp[0] = visIn[joff+0]; dtemp[1] = visIn[joff+1]; 
      dtemp[6] = visIn[joff+3]; dtemp[7] = visIn[joff+4]; 
      /* Special hack for circular feeds with elp/ori solutions
	 zero cross pols and calculate model and correct vis.
	 Don't know why this is necessary */
      if ((me->polType==OBIT_UVPoln_ELORI) && me->circFeed) {
	dtemp[2] = 0.0; dtemp[3] = 0.0;
	dtemp[4] = 0.0; dtemp[5] = 0.0; 
      } else {
	dtemp[2] = visIn[joff+6]; dtemp[3] = visIn[joff+7];
	dtemp[4] = visIn[joff+9]; dtemp[5] = visIn[joff+10];
     }
      /* Now apply calibration by multiplying the inverse Mueller matrix by 
	 the data vector. */
      j = 0;
      if (me->perChan) jndex = (ifoff + choff);
      else             jndex = ifoff;
      /* Is poln cal flagged? */
      if (me->PolCal[jndex]!=fblank) { /* OK */

	/* Now apply calibration by multiplying the inverse Mueller matrix by 
	   the data vector.*/
	MatxVec4Mult(&me->PolCal[jndex], dtemp, ytemp); 
	
	index = 0;
	j = 0;
	/* Reorder - make corrections for cross hand */
	visIn[joff+0]  = ytemp[0];  visIn[joff+1]   = ytemp[1]; 
	visIn[joff+3]  = ytemp[6];  visIn[joff+4]   = ytemp[7]; 
	/* Special hack for circular feeds with elp/ori solutions*/
	if ((me->polType==OBIT_UVPoln_ELORI) && me->circFeed) {
	  visIn[joff+6] += ytemp[2];  visIn[joff+7]  += ytemp[3]; 
	  visIn[joff+9] += ytemp[4];  visIn[joff+10] += ytemp[5];
	} else {
	  visIn[joff+6]  = ytemp[2];  visIn[joff+7]   = ytemp[3]; 
	  visIn[joff+9]  = ytemp[4];  visIn[joff+10]  = ytemp[5];
	}
      } else { /* Bad - flag crosspol */
 	visIn[joff+6] = visIn[joff+7]  = visIn[joff+8]  = 0.0;
 	visIn[joff+9] = visIn[joff+10] = visIn[joff+11] = 0.0;
      }
	
      /* Done if in->numStok < 4) */
      if (in->numStok < 4) {choff += 32; continue;}
      
      /* parallactic angle  - not for circular 'ori-eli' and 'vlbi' */
      if ((me->polType==OBIT_UVPoln_ELORI && me->circFeed) || 
	  (me->polType==OBIT_UVPoln_VLBI)) {
	gr = 1.0;
	gi = 0.0;
      } else if (me->polType==OBIT_UVPoln_ELORI && !me->circFeed) { /* Linear feeds */
	/* Correct for parallactic angle */
	gr = me->curCosPA[ia1] * me->curCosPA[ia2] - me->curSinPA[ia1] * me->curSinPA[ia2];
	gi = me->curCosPA[ia1] * me->curSinPA[ia2] + me->curSinPA[ia1] * me->curCosPA[ia2]; 
      } else { /* others - need parallactic angle correction */
	gr = me->curCosPA[ia1] * me->curCosPA[ia2] - me->curSinPA[ia1] * me->curSinPA[ia2];
	gi = me->curCosPA[ia1] * me->curSinPA[ia2] + me->curSinPA[ia1] * me->curCosPA[ia2];
      }
      
      /* correct RL,LR for parallactic angle and ionospheric Faraday rotation: */
      if (cal!=NULL) {
	loff = (iif - 1) * cal->numLambda + ifreq - 1;
	Lambda2 = cal->Lambda[loff]*cal->Lambda[loff];
	gr1 = gr * cos (Lambda2 * (cal->IFR[ia1] + cal->IFR[ia2])) -  
	  gi * sin (Lambda2 * (cal->IFR[ia1] + cal->IFR[ia2]));
	gi1 = gi * cos (Lambda2 * (cal->IFR[ia1] + cal->IFR[ia2])) +  
	  gr * sin (Lambda2 * (cal->IFR[ia1] + cal->IFR[ia2]));
      } else {
	gr1 = gr;
	gi1 = gi;
      }
      
      /* Correct RL */
      jrl = joff + 2*desc->incs;
      tr = visIn[jrl];
      ti = visIn[jrl+1];
      visIn[jrl]   = tr * gr1 - ti * gi1;
      visIn[jrl+1] = ti * gr1 + tr * gi1;
      
      /* Correct LR */
      jlr = joff + 3*desc->incs;
      tr = visIn[jlr];
      ti = visIn[jlr+1];
      visIn[jlr]   = tr * gr1 + ti * gi1;
      visIn[jlr+1] = ti * gr1 - tr * gi1;
      choff += 32;   /* Channel offset in PolCal */
    } /* end loop over channels L300: */
    ifoff += 32*chdelta;  /* IF offset in PolCal */
  } /* end loop  over IF L400: */
  
} /* end ObitUVCalPolarization */

/**
 * Shutdown Polarization calibration
 * Destroy structures.
 * \param in   Polarization Object.
 * \param err  ObitError stack.
 */
void ObitUVCalPolarizationShutdown (ObitUVCal *in, ObitErr *err)
{
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVCalIsA(in));

  /* delete structure */
  in->polnCal = ObitUVCalPolarizationSUnref(in->polnCal);
} /*  end ObitUVCalPolarizationShutdown */

/**
 * Destroy structure for baseline dependent Calibrate .
 * \param in   Polarization Object.
 * \return NULL
 */
ObitUVCalPolarizationS*
ObitUVCalPolarizationSUnref (ObitUVCalPolarizationS *in)
{
  olong i;

  if (in==NULL) return in;
  
  in->PCal = ObitPolCalListUnref(in->PCal);
  if (in->curPA)    g_free(in->curPA);
  if (in->curCosPA) g_free(in->curCosPA);
  if (in->curSinPA) g_free(in->curSinPA);
  if (in->PolCal)   g_free(in->PolCal);
  if (in->Jones) {
    for (i=0; i<in->numAnt; i++) if (in->Jones[i]) g_free(in->Jones[i]);
    g_free(in->Jones);
  }

  /* basic structure */
   g_free (in);

  return NULL;
} /*  end ObitUVCalPolarizationSUnref */

/*---------------Private functions---------------------------*/

/**
 * Create structure for baseline dependent Calibrate .
 * \param in   Polarization Object.
 * \return newly created object.
 */
static ObitUVCalPolarizationS*
newObitUVCalPolarizationS (ObitUVCal *in)
{
  ObitUVCalPolarizationS* out;

  out = g_malloc0(sizeof(ObitUVCalPolarizationS));

  /* Init pointers */
  out->curPA        = NULL;
  out->curCosPA     = NULL;
  out->curSinPA     = NULL;
  out->PolCal       = NULL;
  out->Jones        = NULL;
  out->PCal         = NULL;

  return out;
} /*  end newObitUVCalPolarizationS */

/**
 * Update ObitUVCalPolarization calibration tables for time time.
 * Calculate a new set of parallactic angles.
 * \param in      Polarization Object.
 * \param UVCal   Basic UV calibration structure
 * \param cal     amp/phase calibration structure with IFR information
 * \param time    desired time in days
 * \param SourID  Source ID, <= 0 -> none.
 * \param SubA    Subarray number
 * \param FreqID  Frequency ID, checked against table
 * \param err     Error Stack.
 */
static void ObitUVCalPolarizationUpdate (ObitUVCalPolarizationS *in, ObitUVCal *UVCal, 
					 ObitUVCalCalibrateS *cal, 
					 ofloat time, olong SourID, olong SubA, olong FreqID, 
					 ObitErr *err)
{
  olong sid, ichan, i;
  odouble Dec, ArrLong, ArrLat, AntLst, HrAng;
  ofloat PA;
  ObitAntennaList *Ant;
  gchar *routine="ObitUVCalPolarizationUpdate";

  /* Check that have calibration */
  if ((UVCal->antennaLists==NULL) || (UVCal->antennaLists[SubA-1]==NULL)) {
    Obit_log_error(err, OBIT_Error, "No polarization Cal info for %s", "Poln Cal");
   }
 
 /* Check FreqID */
  if ((FreqID>0) && (UVCal->antennaLists[SubA-1]->FreqID>0) && 
      (UVCal->antennaLists[SubA-1]->FreqID != FreqID)) {
    Obit_log_error(err, OBIT_Error, "Wrong FreqID %d for Poln Cal for  %s", FreqID, "Poln Cal");
  }
  if (err->error) return; /* Bail out if trouble */

  /* Use central frequency channel, each IF */
  ichan = 0.5 * (in->bChan + in->eChan);

  /* New Source? */
  if (in->curSourID != SourID) {
    sid = MAX (1, SourID) - 1;
    in->curRA     = UVCal->sourceList->SUlist[sid]->RAApp*DG2RAD;
    Dec = UVCal->sourceList->SUlist[sid]->DecApp*DG2RAD;
    in->curCosDec = cos(Dec);
    in->curSinDec = sin(Dec);
  }

  /* Parallactic angles */
  Ant = UVCal->antennaLists[SubA-1];
  /* VLA all the same */
  if (Ant->isVLA ) {
    ArrLong = Ant->ANlist[0]->AntLong;
    ArrLat  = Ant->ANlist[0]->AntLat;
    AntLst = Ant->GSTIAT0 + ArrLong + time*Ant->RotRate;
    HrAng = AntLst - in->curRA;
    PA = atan2 (cos (ArrLat) * sin (HrAng), 
		(sin (ArrLat) * in->curCosDec - cos (ArrLat) * in->curSinDec * cos(HrAng)));
    for (i=0; i<Ant->number; i++) {
      in->curPA[i]    = PA;
      in->curCosPA[i] = cos(PA);
      in->curSinPA[i] = sin(PA);
    }
  } else { /* Not VLA - calculate each */
    for (i=0; i<Ant->number; i++) {
      /* Alt-Az mount and Valid data? */
      if ((Ant->ANlist[i]->AntID>0) && (Ant->ANlist[i]->AntMount==0)) {
	ArrLong = Ant->ANlist[i]->AntLong;
	ArrLat  = Ant->ANlist[i]->AntLat;
	AntLst  = Ant->GSTIAT0 + ArrLong + time*Ant->RotRate;
	HrAng   = AntLst - in->curRA;
	PA = atan2 (cos (ArrLat) * sin (HrAng), 
		    (sin (ArrLat) * in->curCosDec - cos (ArrLat) * in->curSinDec * cos(HrAng)));
	in->curPA[i]    = PA; 
	in->curCosPA[i] = cos(PA);
	in->curSinPA[i] = sin(PA);
      } else { /* Not alt-az or no data */
	in->curPA[i]    = 0.0;
	in->curCosPA[i] = 1.0;
	in->curSinPA[i] = 0.0;
      }
    }
  } /* end non-VLA */

  /* Set inverse Jones matrices if needed */
  /* Do if first time, new subarray or not simple VLA style linear approximation. */
  if ((in->curTime<-100.0) || (in->curSubA  != SubA) || (Ant->polType!=OBIT_UVPoln_Approx)) {
    for (i=0; i<Ant->number; i++) {
      /* Per IF (AIPS AN) or channel (AIPS PD)? */
      if (in->perChan) SetInvJonesCh(in, cal, ichan, i+1, err);  /* Per Channel */
      else             SetInvJonesIF(in, Ant, cal, ichan, i+1);  /* Per IF */
      if (err->error) Obit_traceback_msg (err, routine, in->name);
    }
  }

  /* Save current values */
  in->curTime   = time;
  in->curSourID = SourID;
  in->curSubA   = SubA;
} /* end ObitUVCalPolarizationUpdate */

/**
 * Form IF baseline correction Mueller matrix 
 * R/L Linear model, one matrix per IF only.
 * Use VLBIPol for per channel/IF
 * Use fact that inverse of outer product is outer product of inverse matrices.
 * Matrix assumes data order RR,RL,LR,LL
 * Use method of AIPS/POLSET.FOR
 * \param in      Polarization Object.
 * \param UVCal   Basic UV calibration structure
 * \param SubA    Subarray number (1-rel)
 * \param iant1   First antenna of baseline (1-rel)
 * \param iant2   Second antenna of baseline (1-rel)
 * \param iChan   Central channel to use
 * \param err     Error Stack.
 */
static void LinPolIF(ObitUVCalPolarizationS *in, ObitUVCal *UVCal,olong SubA, 
		     olong iant1, olong iant2, olong iChan, ObitErr *err)
{
  olong i, iif, jif, ia1, ia2, ifoff;
  ofloat D1r[2], D1l[2], D2r[2], D2l[2];
  ObitAntennaList *Ant;
  /* ofloat dbg[4][8];  debug */
  /*gint iii;  debug */

  /* loop thru IFs */
  ifoff = 0; /* Offset in PolCal to beginning of IF matrix */
  ia1 = iant1 - 1;
  ia2 = iant2 - 1;
  Ant = UVCal->antennaLists[SubA-1];
  for (iif=  in->bIF; iif<=in->eIF; iif++) { /* loop 400 */
    jif = iif - 1;
    for (i=0; i<32; i++) in->PolCal[ifoff+i] = 0.0;  /* Zero most terms */
    in->PolCal[ifoff+0]  = 1.0;                      /* diagonal terms */
    in->PolCal[ifoff+10] = 1.0;
    in->PolCal[ifoff+20] = 1.0;
    in->PolCal[ifoff+30] = 1.0;
    
    /* D terms from antenna array */
    D1r[0] =  Ant->ANlist[ia1]->FeedAPCal[jif*2+0];
    D1r[1] =  Ant->ANlist[ia1]->FeedAPCal[jif*2+1];
    D1l[0] =  Ant->ANlist[ia1]->FeedBPCal[jif*2+0];
    D1l[1] =  Ant->ANlist[ia1]->FeedBPCal[jif*2+1];
    D2r[0] =  Ant->ANlist[ia2]->FeedAPCal[jif*2+0];
    D2r[1] =  Ant->ANlist[ia2]->FeedAPCal[jif*2+1];
    D2l[0] =  Ant->ANlist[ia2]->FeedBPCal[jif*2+0];
    D2l[1] =  Ant->ANlist[ia2]->FeedBPCal[jif*2+1];
    in->PolCal[ifoff+2] = -0.5 * (D1r[0] + D2l[0]);
    in->PolCal[ifoff+3] = -0.5 * (D1r[1] - D2l[1]);
    in->PolCal[ifoff+3*8+2] = in->PolCal[ifoff+2];
    in->PolCal[ifoff+3*8+3] = in->PolCal[ifoff+3];
    
    in->PolCal[ifoff+4] = -0.5 * (D1l[0] + D2r[0]);
    in->PolCal[ifoff+5] = -0.5 * (D1l[1] - D2r[1]);
    in->PolCal[ifoff+3*8+4] = in->PolCal[ifoff+4];
    in->PolCal[ifoff+3*8+5] = in->PolCal[ifoff+5];
    ifoff += 32; /* Offset in PolCal to beginning of IF matrix */
  } /* end loop over IF */

} /* end  LinPolIF */

/**
 * Form channel/IF baseline correction Mueller matrix 
 * R/L Linear model, one matrix per channel/if.
 * Use VLBIPol for per channel/IF
 * Use fact that inverse of outer product is outer product of inverse matrices.
 * Matrix assumes data order RR,RL,LR,LL
 * Use method of AIPS/POLSET.FOR
 * \param in      Polarization Object.
 * \param UVCal   Basic UV calibration structure
 * \param iant1   First antenna of baseline (1-rel)
 * \param iant2   Second antenna of baseline (1-rel)
 * \param iChan   Central channel to use
 * \param err     Error Stack.
 */
static void LinPolCh(ObitUVCalPolarizationS *in, ObitUVCal *UVCal,
		     olong iant1, olong iant2, olong iChan, ObitErr *err)
{
  olong i, iif, ia1, ia2, ich, nch, jndex, kndx;
  ofloat D1r[2], D1l[2], D2r[2], D2l[2];
  ObitPolCalList *PCal;
  /* ofloat dbg[4][8];  debug */
  /*gint iii;  debug */

  /* loop thru IFs */
  jndex = 0; /* Offset in PolCal to beginning of IF matrix */
  ia1   = iant1 - 1;
  ia2   = iant2 - 1;
  nch   = in->eChan - in->bChan + 1;
  PCal  = in->PCal;
  for (iif= in->bIF; iif<=in->eIF; iif++) { /* loop 400 */
    /* Loop over channel */
    for (ich=in->bChan; ich<=in->eChan; ich++) {
      kndx = 4*((iif-1)*nch + (ich-1));
      
      for (i=0; i<32; i++) in->PolCal[jndex+i] = 0.0;  /* Zero most terms */
      in->PolCal[jndex+0]  = 1.0;                      /* diagonal terms */
      in->PolCal[jndex+10] = 1.0;
      in->PolCal[jndex+20] = 1.0;
      in->PolCal[jndex+30] = 1.0;
      
      /* D terms from antenna array */
      D1r[0] =  PCal->ANlist[ia1][kndx+0];
      D1r[1] =  PCal->ANlist[ia1][kndx+1];
      D1l[0] =  PCal->ANlist[ia1][kndx+2];
      D1l[1] =  PCal->ANlist[ia1][kndx+3];
      D2r[0] =  PCal->ANlist[ia2][kndx+0];
      D2r[1] =  PCal->ANlist[ia2][kndx+1];
      D2l[0] =  PCal->ANlist[ia2][kndx+2];
      D2l[1] =  PCal->ANlist[ia2][kndx+3];
      in->PolCal[jndex+2] = -0.5 * (D1r[0] + D2l[0]);
      in->PolCal[jndex+3] = -0.5 * (D1r[1] - D2l[1]);
      in->PolCal[jndex+3*8+2] = in->PolCal[jndex+2];
      in->PolCal[jndex+3*8+3] = in->PolCal[jndex+3];
      
      in->PolCal[jndex+4] = -0.5 * (D1l[0] + D2r[0]);
      in->PolCal[jndex+5] = -0.5 * (D1l[1] - D2r[1]);
      in->PolCal[jndex+3*8+4] = in->PolCal[jndex+4];
      in->PolCal[jndex+3*8+5] = in->PolCal[jndex+5];
      jndex += 32; /* Offset in PolCal to beginning of IF matrix */
    } /* end loop over channel */
  } /* end loop over IF */

} /* end  LinPolCh */

/**
 * Form baseline inverse Mueller matrix from antenna inverse Jones matrices
 * Use fact that inverse of outer product is outer product of inverse matrices.
  * Elipticity/Orientation model, one matrix per IF or channel/IF.
 * Matrix assumes data order RR,RL,LR,LL
 * \param in      Polarization Object.
 * \param iant1   First antenna of baseline
 * \param iant2   Second antenna of baseline
 * \param iChan   Central channel to use
 * \param err     Error Stack.
 */
static void OriPol(ObitUVCalPolarizationS *in, olong iant1, olong iant2, olong iChan, ObitErr *err)
{
  olong ia1, ia2, iif, ich, bChan, eChan, ioff, jndex;

  /* Per Channels? */
  if (in->perChan) {bChan = in->bChan; eChan = in->eChan; }
  else             {bChan = 1;         eChan = bChan; }

  /* loop thru IFs */
  ioff  = 0; /* Offset in Jones to beginning of IF/channel matrix */
  jndex = 0; /* Offset in PolCal to beginning of IF/channel matrix */
  ia1 = iant1 - 1;
  ia2 = iant2 - 1;
  for (iif=  in->bIF; iif<=in->eIF; iif++) { /* loop 400 */
 
    /* Loop over channel */
    for (ich=bChan; ich<=eChan; ich++) {
      
      /* Muller matrix = outer product iant1 * conjg(iant2), for order RR,RL,LR,LL */
      MatxOuter (&in->Jones[ia1][ioff], &in->Jones[ia2][ioff], &in->PolCal[jndex]);

      jndex += 32; /* Offset in PolCal to beginning of channel/IF matrix */
      ioff  +=  8; /* Offset in Jones to beginning of  channel/IF  */
    } /* end channel loop */
  } /* end loop over IF */
} /* end OriPol  */

/**
 * Form baseline inverse Mueller matrix from antenna inverse Jones matrices
 * R/L Linear model resolved source, one matrix per IF or channel/IF.
 * Use fact that inverse of outer product is outer product of inverse matrices.
 * Matrix assumes data order RR,RL,LR,LL
 * \param in      Polarization Object.
 * \param iant1   First antenna of baseline
 * \param iant2   Second antenna of baseline
 * \param iChan   Central channel to use
 * \param err     Error Stack.
 */
static void VLBIPol(ObitUVCalPolarizationS *in, olong iant1, olong iant2, olong iChan, ObitErr *err)
{
  olong iif, ia1, ia2, ioff, ich, bChan, eChan, jndex;

  /* Per Channels? */
  if (in->perChan) {bChan = in->bChan; eChan = in->eChan; }
  else             {bChan = 1;         eChan = bChan; }

  /* loop thru IFs */
  jndex = 0; /* Offset in PolCal to beginning of IF matrix */
  ioff  = 0; /* Offset in Jones to beginning of IF matrix */
  ia1 = iant1 - 1;
  ia2 = iant2 - 1;
  for (iif=  in->bIF; iif<=in->eIF; iif++) { /* loop 400 */

    /* Loop over channel */
    for (ich=bChan; ich<=eChan; ich++) {
      
      /* outer product iant1 * conjg(iant2), for order RR,RL,LR,LL */
      MatxOuter (&in->Jones[ia1][ioff], &in->Jones[ia2][ioff], &in->PolCal[jndex]);
      
      ioff  +=  8; /* Offset in Jones to beginning of IF matrix */
      jndex += 32; /* Offset in PolCal to beginning of IF matrix */
    } /* end channel loop */
  } /* end loop over IF */
} /* end VLBIPol  */

/**
 * Form baseline inverse Mueller matrix from antenna inverse Jones matrices
 * X/Y Linear model, one matrix per IF or channel/IF.
 * Resulting matrix must also include transformations to convert to RR,LL,RL,LR.
 * NYI.
 * \param in      Polarization Object.
 * \param iant1   First antenna of baseline
 * \param iant2   Second antenna of baseline
 * \param iChan   Central channel to use
 * \param err     Error Stack.
 */
static void LXYPol(ObitUVCalPolarizationS *in, olong iant1, olong iant2, olong iChan, ObitErr *err)
{
  g_error ("LXYPol: X/Y polarization not yet implemented ");
} /* end LXYPol */

/**
 * Form Antenna inverse Jones matrices per IF from AIPS AN table
 * For the OBIT_UVPoln_Approx model, the calibration does not need to include
 * parallactic angle (and IFR)
 * Use fact that inverse of outer product is outer product of inverse matrices.
 * \param in      Polarization Object.
 * \param Ant     Antenna list object
 * \param cal     Amp/phase calibration object for IFR calibration.
 * \param iant    Antenna number
 */
static void SetInvJonesIF(ObitUVCalPolarizationS *in, ObitAntennaList *Ant, 
			  ObitUVCalCalibrateS *cal, olong iChan, olong iant)
{
  ofloat Dr[2]={0.0,0.0}, Dl[2]={0.0,0.0}, d;
  ofloat Det[2] = {1.0, 0.0};
  ofloat Jones[8] = {1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0};
  ofloat elp_r, elp_l, ori_r, ori_l, angle[6], sina[6], cosa[6];
  ofloat root2, fblank = ObitMagicF();
  ofloat rotate=0.0, chi=0.0, crot, srot, PD, temp[8];
  olong iif, jndx, loff, refAnt, i, SubA;
  gboolean doJones=TRUE;
  ocomplex RS, RD, LS, LD, PA, PAc, PRref, PLref, ct;
  ObitPolCalList *PCal = in->PCal;

  SubA = MAX (1, in->curSubA);
  refAnt = Ant->polRefAnt;
  root2  = 1.0 / sqrt(2.0);

  /* WHERE IS THIS FROM??? */
  PD = 0.0;

  /* parallactic angle */
  COMPLEX_EXP (PA,  2*in->curPA[iant-1]);
  COMPLEX_CONJUGATE (PAc, PA);

  /* Loop over IFs (index 0 rel) */
  jndx = 0;
  for (iif= in->bIF-1; iif<=in->eIF-1; iif++) {

    /* Set Jones matrix by type. */
    /* Each entry is a 2x2 complex matrix - it is way beyond c to represent this */
    switch (Ant->polType) {
    case OBIT_UVPoln_ELORI:   /* Elipticity-orientation - SPECIAL CASE */
      /* In AIPS AN table, these are stored as (elipticity,Orientation) in 
	 POLCA and POLCB */
      doJones = FALSE;  /* Jones matrix terms computed here */
      elp_r = Ant->ANlist[iant-1]->FeedAPCal[2*iif+0];
      if (refAnt>0) ori_r = Ant->ANlist[refAnt-1]->FeedAPCal[2*iif+1];
      else          ori_r = 0.0;
      elp_l = Ant->ANlist[iant-1]->FeedBPCal[2*iif+0];
      if (refAnt>0) ori_l = Ant->ANlist[refAnt-1]->FeedBPCal[2*iif+1];
      else          ori_l = 0.0;
      
      if ((elp_r!=fblank) && (elp_l!=fblank)) {
	/* OK */
	/* for circular feeds */
	if (in->circFeed) {
	  
	  /* Sines/cosines */
	  angle[0] = elp_r; angle[1] = 2*ori_r; angle[2] = elp_l; angle[3] = -2*ori_l;
	  if (PCal->polRefAnt>0) angle[4] = Ant->ANlist[PCal->polRefAnt-1]->FeedAPCal[2*iif+1];
	  else                   angle[4] = 0.0;
	  /* RLPhaseDiff should have been applied in BP table 
	     if (PCal->polRefAnt>0) angle[5] = -Ant->ANlist[PCal->polRefAnt-1]->FeedBPCal[2*iif+1] + Ant->RLPhaseDiff[iif];
	     else                   angle[5] =  Ant->RLPhaseDiff[iif];*/
	  if (PCal->polRefAnt>0) angle[5] = -Ant->ANlist[PCal->polRefAnt-1]->FeedBPCal[2*iif+1];
	  else                   angle[5] =  0.0;
	  ObitSinCosVec(6, angle, sina, cosa);
	  
	  /* Complex terms */
	  COMPLEX_SET (RS, root2*(cosa[0] + sina[0]), 0.0);
	  COMPLEX_SET (RD, root2*(cosa[0] - sina[0]) * cosa[1], root2*(cosa[0] - sina[0]) * sina[1]);
	  COMPLEX_SET (LS, root2*(cosa[2] + sina[2]) * cosa[3], root2*(cosa[2] + sina[2]) * sina[3]);
	  COMPLEX_SET (LD, root2*(cosa[2] - sina[2]), 0.0);
	  COMPLEX_SET (PRref, cosa[4], sina[4]);
	  COMPLEX_SET (PLref, cosa[5], sina[5]);
	  
	  /* Jones matrix: 
	     RS * PRref         RD * PA * PRref
	     LS * PAc * PLref   LD * PLref   
	  */
	  COMPLEX_MUL2 (ct, RS, PRref);
	  Jones[0] = ct.real;
	  Jones[1] = ct.imag;
	  COMPLEX_MUL3 (ct, RD, PA, PRref);
	  Jones[2] = ct.real;
	  Jones[3] = ct.imag;
	  COMPLEX_MUL3 (ct, LS, PAc, PLref);
	  Jones[4] = ct.real;
	  Jones[5] = ct.imag;
	  COMPLEX_MUL2 (ct, LD, PLref);
	  Jones[6] = ct.real;
	  Jones[7] = ct.imag;
	  
	} else {  /* Linear */
	  /* Jones matrix: 
	     i cos(pi/4+elp_r)exp(-i ori_r)exp(i -chi)  sin(pi/4+elp_r)exp(+i ori_r)exp(i chi)
	    -i sin(pi/4-elp_l)exp(-i ori_l)exp(i -chi) -cos(pi/4-elp_l)exp(+i ori_l)exp(i chi)
	  */
	  /* Sines/cosines */
	  chi = in->curPA[iant-1];  /* Parallactic angle */
	  angle[0] = G_PI*0.25+elp_r; angle[1] = G_PI*0.25-elp_l;
	  angle[2] = ori_r+chi;       angle[3] = ori_l+chi;
	  ObitSinCosVec(4, angle, sina, cosa);
	  Jones[0] =  cosa[0]*sina[2];
	  Jones[1] =  cosa[0]*cosa[2];
	  Jones[2] =  sina[0]*cosa[2];
	  Jones[3] =  sina[0]*sina[2];
	  Jones[4] = -sina[1]*sina[3];  /* Flip sign of "Y" */
	  Jones[5] = -sina[1]*cosa[3];
	  Jones[6] = -cosa[1]*cosa[3];
	  Jones[7] = -cosa[1]*sina[3];
	} /* end linear feeds */
	/* Also need (inverse of) determinant */
	Det[0] = (Jones[0]*Jones[6] - Jones[1]*Jones[7]) - (Jones[2]*Jones[4] - Jones[3]*Jones[5]);
	Det[1] = (Jones[0]*Jones[7] + Jones[1]*Jones[6]) - (Jones[2]*Jones[5] + Jones[3]*Jones[4]);
	/* Inverse of determinant */
	d = Det[0]*Det[0] + Det[1]*Det[1];
	if (d!=0.0) d = 1.0 / d;
	else d = 1.0;
	Det[0] *=  d;
	Det[1] *= -d;
      } else { /* bad */
	for (i=0; i<8; i++) Jones[i] = fblank;
	Det[0] = 1.0; Det[1] = 0.0;
	rotate = 0.0;
     }

      break;

    case OBIT_UVPoln_Approx:  /* R/L Linear D-term approximation */
      Dr[0] =  Ant->ANlist[iant-1]->FeedAPCal[2*iif+0];
      Dr[1] =  Ant->ANlist[iant-1]->FeedAPCal[2*iif+1];
      Dl[0] =  Ant->ANlist[iant-1]->FeedBPCal[2*iif+0];
      Dl[1] =  Ant->ANlist[iant-1]->FeedBPCal[2*iif+1];
     /* Don't need rotation by parallactic angle, IFR */
      rotate = 0.0;
      break;

    case OBIT_UVPoln_VLBI:    /* R/L Linear D-term approximation for resolved sources */
      Dr[0] =  Ant->ANlist[iant-1]->FeedAPCal[2*iif+0];
      Dr[1] =  Ant->ANlist[iant-1]->FeedAPCal[2*iif+1];
      Dl[0] =  Ant->ANlist[iant-1]->FeedBPCal[2*iif+0];
      Dl[1] =  Ant->ANlist[iant-1]->FeedBPCal[2*iif+1];
      /* debug zero PA
      in->curPA[iant-1] = 0.0;
      in->curCosPA[iant-1] = 1.0;
      in->curSinPA[iant-1] = 0.0; */

      rotate = in->curPA[iant-1];
      if ((cal!=NULL) && (cal->IFR[iant-1] != fblank)) {
	loff = (iif - 1) * cal->numLambda + iChan - 1;
	rotate += cal->Lambda[loff]*cal->Lambda[loff] * cal->IFR[iant-1];
      }

      /* Rotate D terms in the way the Data were rotated */
      rotate *= 2.0;
      crot = cos(rotate);
      srot = sin(rotate);
      temp[0] = Dr[0]; temp[1] = Dr[1];
      Dr[0] = temp[0] * crot - temp[1] * srot;
      Dr[1] = temp[0] * srot + temp[1] * crot;
      temp[0] = Dl[0]; temp[1] = Dl[1];
      Dl[0] =  temp[0] * crot + temp[1] * srot;
      Dl[1] = -temp[0] * srot + temp[1] * crot;
      rotate = 0.0;
      break;

    case OBIT_UVPoln_XYLin:   /* X/Y Linear D-term approximation  */
      Dr[0] =  Ant->ANlist[iant-1]->FeedAPCal[2*iif+0];
      Dr[1] =  Ant->ANlist[iant-1]->FeedAPCal[2*iif+1];
      Dl[0] =  Ant->ANlist[iant-1]->FeedBPCal[2*iif+0];
      Dl[1] =  Ant->ANlist[iant-1]->FeedBPCal[2*iif+1];
      break;
    default:
      g_assert_not_reached(); /* unknown, barf */
    }; /* end switch */
    
    if (doJones) {
      /* Calculate Jones Matrix */
      Jones[0] = 1.0;
      Jones[1] = 0.0;
      Jones[2] = Dr[0];
      Jones[3] = Dr[1];
      Jones[4] = Dl[0];
      Jones[5] = Dl[1];
      Jones[6] = 1.0;
      Jones[7] = 0.0;
      
      /* Rotate elements by parallactic angle/IFR if needed */
      if (fabs(rotate) > 0.001) {
	for (i=0; i<8; i++) temp[i] = Jones[i]; /* copy of Jones matrix */
	crot = cos(rotate);
	srot = sin(rotate);
	Jones[0] = temp[0] * crot - temp[1] * srot;
	Jones[1] = temp[0] * srot + temp[1] * crot;
	Jones[2] = temp[2] * crot - temp[3] * srot;
	Jones[3] = temp[2] * srot + temp[3] * crot;
	Jones[4] = temp[4] * crot - temp[5] * srot;
	Jones[5] = temp[4] * srot + temp[5] * crot;
	Jones[6] = temp[6] * crot - temp[7] * srot;
	Jones[7] = temp[6] * srot + temp[7] * crot;
	
	/* Also need (inverse of) determinant */
	Det[0] = (Jones[0]*Jones[6] - Jones[1]*Jones[7]) - (Jones[2]*Jones[4] - Jones[3]*Jones[5]);
	Det[1] = (Jones[0]*Jones[7] + Jones[1]*Jones[6]) - (Jones[2]*Jones[5] + Jones[3]*Jones[4]);
	/* Inverse of determinant */
	d = Det[0]*Det[0] + Det[1]*Det[1];
	if (d!=0.0) d = 1.0 / d;
	else d = 1.0;
	Det[0] *=  d;
	Det[1] *= -d;
	
      } else { /* only need inverse of determinant */
	Det[0] =  1.0 - Dr[0]*Dl[0] + Dr[1]*Dl[1];
	Det[1] =  -Dr[0]*Dl[1] - Dr[1]*Dl[0];
	/* Inverse of determinant */
	d = Det[0]*Det[0] + Det[1]*Det[1];
	if (d!=0.0) d = 1.0 / d;
	else d = 1.0;
	Det[0] *=  d;
	Det[1] *= -d;
      }
    } /* end doJones */
    
    /* invert matrix if valid */
    if (Jones[0]!=fblank) {
      in->Jones[iant-1][jndx+6] =   Jones[0] * Det[0] - Jones[1] * Det[1];
      in->Jones[iant-1][jndx+7] =   Jones[0] * Det[1] + Jones[1] * Det[0];
      in->Jones[iant-1][jndx+2] = -(Jones[2] * Det[0] - Jones[3] * Det[1]);
      in->Jones[iant-1][jndx+3] = -(Jones[2] * Det[1] + Jones[3] * Det[0]);
      in->Jones[iant-1][jndx+4] = -(Jones[4] * Det[0] - Jones[5] * Det[1]);
      in->Jones[iant-1][jndx+5] = -(Jones[4] * Det[1] + Jones[5] * Det[0]);
      in->Jones[iant-1][jndx+0] =   Jones[6] * Det[0] - Jones[7] * Det[1];
      in->Jones[iant-1][jndx+1] =   Jones[6] * Det[1] + Jones[7] * Det[0];
    } else { /* bad */
      for (i=0; i<8; i++) in->Jones[iant-1][jndx+i] = fblank;
    }
    
    jndx += 8; /* Increment index for next IF */
  } /* end loop over IFs */
} /* end SetInvJonesIF */

/**
 * Form Antenna inverse Jones matrices per channel/IF from AIPS PD table
 * For the OBIT_UVPoln_Approx model, the calibration does not need to include
 * parallactic angle (and IFR)
 * Use fact that inverse of outer product is outer product of inverse matrices.
 * \param in      Polarization Object.
 * \param PCal    Poln Cal list object
 * \param cal     Amp/phase calibration object for IFR calibration.
 * \param iant    Antenna number 1-rel
 */
static void SetInvJonesCh(ObitUVCalPolarizationS *in, ObitUVCalCalibrateS *cal, 
			  olong iChan, olong iant, ObitErr *err)
{
  ofloat Jones[8], Dr[2]={0.0,0.0}, Dl[2]={0.0,0.0}, Det[2], d;
  ofloat elp_r, elp_l, ori_r, ori_l, PD=0.0, angle[6], sina[6], cosa[6];
  ofloat root2, fblank = ObitMagicF();
  ofloat rotate=0.0, chi=0.0, crot, srot, temp[8];
  olong iif, jndx, loff, refAnt, i, SubA, ich, nch, nchAll, kndx, lndx, ia;
  gboolean doJones=TRUE;
  ObitPolCalList *PCal = in->PCal;
  ocomplex RS, RD, LS, LD, PA, PAc, PRref, PLref, ct;
  /*ocomplex CX, SX, CY, SY;*/
  gchar *routine="SetInvJonesCh";

  if (err->error) return;  /* prior error? */
  g_assert (PCal!=NULL);

  /* Check number of channels */
  Obit_return_if_fail((in->numChan==PCal->numChan), err, 
		      "%s: Unequal channels %d != %d for %s", 
		      routine, in->numChan, PCal->numChan, in->name);

  /* Check number of IFs */
  Obit_return_if_fail((in->numIF==PCal->numIF), err, 
		      "%s: Unequal IFs %d != %d for %s", 
		      routine, in->numIF, PCal->numIF, in->name);

  SubA   = MAX (1, in->curSubA);
  refAnt = PCal->polRefAnt;
  nchAll = in->numChan;                   /* Number of channels in input */
  nch    = in->eChan - in->bChan + 1;     /* Number of channels selected */
  ia     = iant - 1;
  root2  = 1.0 / sqrt(2.0);

  /* parallactic angle */
  COMPLEX_EXP (PA, 2*in->curPA[ia]);
  COMPLEX_CONJUGATE (PAc, PA);

  /* Loop over IFs (index 0 rel) */
  jndx = 0;
  for (iif=in->bIF-1; iif<=in->eIF-1; iif++) {

    /* Loop over channel */
    for (ich=0; ich<nch; ich++) {
      kndx = 4*(iif*nchAll + ich + (in->bChan-1));
      lndx = iif*nchAll + ich + (in->bChan-1);
      
      /* Set Jones matrix by type. */
      doJones = TRUE; 
      /* Each entry is a 2x2 complex matrix - it is way beyond c to represent this */
      switch (PCal->polType) {
      case OBIT_UVPoln_ELORI:   /* Elipticity-orientation  */
	/* In AIPS AN table, these are stored as (ellipticity,Orientation) in 
	   POLCA and POLCB in AN table or in PD table as Elip_R/X, Ori_R/X, Elip_L/Y, Ori_L/Y */
	doJones = FALSE;  /* Jones matrix terms computed here */
	elp_r = PCal->ANlist[ia][kndx+0];
	ori_r = PCal->ANlist[ia][kndx+1];
	elp_l = PCal->ANlist[ia][kndx+2];
	ori_l = PCal->ANlist[ia][kndx+3];
	PD    = PCal->RLPhaseDiff[lndx]*DG2RAD;  /* R-L phase difference */

	if ((elp_r!=fblank) && (elp_l!=fblank)) {
	  /* OK */
	  /* for circular feeds */
	  if (in->circFeed) {
	    
	    /* Sines/cosines */
	    angle[0] = elp_r; angle[1] = 2*ori_r; angle[2] = elp_l; angle[3] = -2*ori_l;
	    if (PCal->polRefAnt>0) angle[4] =  PCal->ANlist[PCal->polRefAnt-1][kndx+1];
	    else                   angle[4] = 0.0;
	    if (PCal->polRefAnt>0) angle[5] = -PCal->ANlist[PCal->polRefAnt-1][kndx+3] + PD;
	    else                   angle[5] = PD;
	    ObitSinCosVec(6, angle, sina, cosa);
	    
	    /* Complex terms */
	    COMPLEX_SET (RS, root2*(cosa[0] + sina[0]), 0.0);
	    COMPLEX_SET (RD, root2*(cosa[0] - sina[0]) * cosa[1], root2*(cosa[0] - sina[0]) * sina[1]);
	    COMPLEX_SET (LS, root2*(cosa[2] + sina[2]) * cosa[3], root2*(cosa[2] + sina[2]) * sina[3]);
	    COMPLEX_SET (LD, root2*(cosa[2] - sina[2]), 0.0);
	    COMPLEX_SET (PRref, cosa[4], sina[4]);
	    COMPLEX_SET (PLref, cosa[5], sina[5]);
	    
	    /* Jones matrix: 
	       RS * PRref         RD * PA * PRref
	       LS * PAc * PLref   LD * PLref   
	    */
	    COMPLEX_MUL2 (ct, RS, PRref);
	    Jones[0] = ct.real;
	    Jones[1] = ct.imag;
	    COMPLEX_MUL3 (ct, RD, PA, PRref);
	    Jones[2] = ct.real;
	    Jones[3] = ct.imag;
	    COMPLEX_MUL3 (ct, LS, PAc, PLref);
	    Jones[4] = ct.real;
	    Jones[5] = ct.imag;
	    COMPLEX_MUL2 (ct, LD, PLref);
	    Jones[6] = ct.real;
	    Jones[7] = ct.imag;
	  } else {  /* Linear */
	    /* Jones matrix: 
	       i cos(pi/4+elp_r)exp(-i ori_r)exp(i -chi)  sin(pi/4+elp_r)exp(+i ori_r)exp(i chi)
	       i sin(pi/4-elp_l)exp(-i ori_l)exp(i -chi)  cos(pi/4-elp_l)exp(+i ori_l)exp(i chi)
	    */
	    /* Sines/cosines */
	    chi = in->curPA[ia];  /* Parallactic angle */
	    /* angle[0] = G_PI*0.25+elp_r; angle[1] = G_PI*0.25-elp_l;
	       angle[2] = ori_r+chi;       angle[3] = ori_l+chi;
	       Obit SinCosVec(4, angle, sina, cosa);
	       Jones[0] =  cosa[0]*sina[2]; Jones[1] =  cosa[0]*cosa[2];
	       Jones[2] =  sina[0]*cosa[2]; Jones[3] =  sina[0]*sina[2];
	       Jones[4] =  sina[1]*sina[3]; Jones[5] =  sina[1]*cosa[3];
	       Jones[6] =  cosa[1]*cosa[3]; Jones[7] =  cosa[1]*sina[3];*/
	    /* Jones matrix: 
	       cos(pi/4+elp_r)exp(-j ori_r)  sin(pi/4+elp_r)exp(+i ori_r)
	       sin(pi/4-elp_l)exp(-j ori_l)  cos(pi/4-elp_l)exp(+i ori_l) 
	       Not a function of chi so only calculate once */
	    if (in->Jones[ia][jndx]!=0.0) return;  /* Already done? */
	    angle[0] = G_PI*0.25+elp_r; angle[1] = G_PI*0.25-elp_l;
	    angle[2] = ori_r;           angle[3] = ori_l;
	    ObitSinCosVec(4, angle, sina, cosa);
	    Jones[0] =  cosa[0]*cosa[2]; Jones[1] = -cosa[0]*sina[2];
	    Jones[2] =  sina[0]*cosa[2]; Jones[3] =  sina[0]*sina[2];
	    Jones[4] =  sina[1]*cosa[3]; Jones[5] = -sina[1]*sina[3];
	    Jones[6] =  cosa[1]*cosa[3]; Jones[7] =  cosa[1]*sina[3];
	  } /* end linear feeds */
	  
	  /* Also need (inverse of) determinant */
	  Det[0] = (Jones[0]*Jones[6] - Jones[1]*Jones[7]) - (Jones[2]*Jones[4] - Jones[3]*Jones[5]);
	  Det[1] = (Jones[0]*Jones[7] + Jones[1]*Jones[6]) - (Jones[2]*Jones[5] + Jones[3]*Jones[4]);
	  /* Inverse of determinant */
	  d = Det[0]*Det[0] + Det[1]*Det[1];
	  if (d!=0.0) d = 1.0 / d;
	  else d = 1.0;
	  Det[0] *=  d;
	  Det[1] *= -d;
	} else { /* bad */
	  for (i=0; i<8; i++) Jones[i] = fblank;
	  Det[0] = 1.0; Det[1] = 0.0;
	  rotate = 0.0;
	}
	break;
	
      case OBIT_UVPoln_Approx:  /* R/L Linear D-term approximation */
	Dr[0] =  PCal->ANlist[ia][kndx+0];
	Dr[1] =  PCal->ANlist[ia][kndx+1];
	Dl[0] =  PCal->ANlist[ia][kndx+2];
	Dl[1] =  PCal->ANlist[ia][kndx+3];
	/* Don't need rotation by parallactic angle, IFR */
	rotate = 0.0;
	break;
	
      case OBIT_UVPoln_VLBI:    /* R/L Linear D-term approximation for resolved sources */
	Dr[0] =  PCal->ANlist[ia][kndx+0];
	Dr[1] =  PCal->ANlist[ia][kndx+1];
	Dl[0] =  PCal->ANlist[ia][kndx+2];
	Dl[1] =  PCal->ANlist[ia][kndx+3];
	/* debug zero PA
	   in->curPA[ia] = 0.0;
	   in->curCosPA[ia] = 1.0;
	   in->curSinPA[ia] = 0.0; */
	
	rotate = in->curPA[ia];
	if ((cal!=NULL) && (cal->IFR[ia] != fblank)) {
	  loff = (iif - 1) * cal->numLambda + iChan - 1;
	  rotate += cal->Lambda[loff]*cal->Lambda[loff] * cal->IFR[ia];
	}
	
	/* Rotate D terms in the way the Data were rotated */
	rotate *= 2.0;
	crot = cos(rotate);
	srot = sin(rotate);
	temp[0] = Dr[0]; temp[1] = Dr[1];
	Dr[0] = temp[0] * crot - temp[1] * srot;
	Dr[1] = temp[0] * srot + temp[1] * crot;
	temp[0] = Dl[0]; temp[1] = Dl[1];
	Dl[0] =  temp[0] * crot + temp[1] * srot;
	Dl[1] = -temp[0] * srot + temp[1] * crot;
	rotate = 0.0;
	break;
	
      case OBIT_UVPoln_XYLin:   /* X/Y Linear D-term approximation  */
	Dr[0] =  PCal->ANlist[ia][kndx+0];
	Dr[1] =  PCal->ANlist[ia][kndx+1];
	Dl[0] =  PCal->ANlist[ia][kndx+2];
	Dl[1] =  PCal->ANlist[ia][kndx+3];
	break;
      default:
	g_assert_not_reached(); /* unknown, barf */
      }; /* end switch */
      
      /* Need to compute Jones matrix and inverse determinant? */
      if (doJones) {
	/* Calculate Jones Matrix */
	Jones[0] = 1.0;
	Jones[1] = 0.0;
	Jones[2] = Dr[0];
	Jones[3] = Dr[1];
	Jones[4] = Dl[0];
	Jones[5] = Dl[1];
	Jones[6] = 1.0;
	Jones[7] = 0.0;
	
	/* Rotate elements by parallactic angle/IFR if needed */
	if (fabs(rotate) > 0.001) {
	  for (i=0; i<8; i++) temp[i] = Jones[i]; /* copy of Jones matrix */
	  crot = cos(rotate);
	  srot = sin(rotate);
	  Jones[0] = temp[0] * crot - temp[1] * srot;
	  Jones[1] = temp[0] * srot + temp[1] * crot;
	  Jones[2] = temp[2] * crot - temp[3] * srot;
	  Jones[3] = temp[2] * srot + temp[3] * crot;
	  Jones[4] = temp[4] * crot - temp[5] * srot;
	  Jones[5] = temp[4] * srot + temp[5] * crot;
	  Jones[6] = temp[6] * crot - temp[7] * srot;
	  Jones[7] = temp[6] * srot + temp[7] * crot;
	  
	  /* Also need (inverse of) determinant */
	  Det[0] = (Jones[0]*Jones[6] - Jones[1]*Jones[7]) - (Jones[2]*Jones[4] - Jones[3]*Jones[5]);
	  Det[1] = (Jones[0]*Jones[7] + Jones[1]*Jones[6]) - (Jones[2]*Jones[5] + Jones[3]*Jones[4]);
	  /* Inverse of determinant */
	  d = Det[0]*Det[0] + Det[1]*Det[1];
	  if (d!=0.0) d = 1.0 / d;
	  else d = 1.0;
	  Det[0] *=  d;
	  Det[1] *= -d;
	  
	} else { /* only need inverse of determinant */
	  Det[0] =  1.0 - Dr[0]*Dl[0] + Dr[1]*Dl[1];
	  Det[1] =  -Dr[0]*Dl[1] - Dr[1]*Dl[0];
	  /* Inverse of determinant */
	  d = Det[0]*Det[0] + Det[1]*Det[1];
	  if (d!=0.0) d = 1.0 / d;
	  else d = 1.0;
	  Det[0] *=  d;
	  Det[1] *= -d;
	}
	
     } /* end doJones */
      /* invert matrix if valid */
      if (Jones[0]!=fblank) {
	in->Jones[ia][jndx+6] =   Jones[0] * Det[0] - Jones[1] * Det[1];
	in->Jones[ia][jndx+7] =   Jones[0] * Det[1] + Jones[1] * Det[0];
	in->Jones[ia][jndx+2] = -(Jones[2] * Det[0] - Jones[3] * Det[1]);
	in->Jones[ia][jndx+3] = -(Jones[2] * Det[1] + Jones[3] * Det[0]);
	in->Jones[ia][jndx+4] = -(Jones[4] * Det[0] - Jones[5] * Det[1]);
	in->Jones[ia][jndx+5] = -(Jones[4] * Det[1] + Jones[5] * Det[0]);
	in->Jones[ia][jndx+0] =   Jones[6] * Det[0] - Jones[7] * Det[1];
	in->Jones[ia][jndx+1] =   Jones[6] * Det[1] + Jones[7] * Det[0];
      } else { /* bad */
        for (i=0; i<8; i++) in->Jones[ia][jndx+i] = fblank;
      }
      jndx += 8; /* Increment index for next channel/IF */
   } /* end channel loop */
  } /* end loop over IFs */
} /* end SetInvJonesCh */

  /** Private: 4x4 complex matrix * 4x1 complex vector multiply */
  static void MatxVec4Mult(ofloat* in1, ofloat* in2, ofloat* out)
{
  olong ic, ii, n=4;
  ofloat sumr, sumi;
  for (ic=0; ic<4; ic++) {
    sumr = sumi = 0.0;
    for (ii=0; ii<4; ii++) {
      sumr += in1[(ic*n+ii)*2]*in2[ii*2]   - in1[(ic*n+ii)*2+1]*in2[ii*2+1];
      sumi += in1[(ic*n+ii)*2]*in2[ii*2+1] + in1[(ic*n+ii)*2+1]*in2[ii*2];
    }
    out[ic*2]   = sumr;
    out[ic*2+1] = sumi;
  }
} /* end MatxVec4Mult */

/** Private: Muller matrix from outer product of Jones matrices 
   Supports blanking */
static void MatxOuter(ofloat* in1, ofloat* in2, ofloat* out)
{
  olong i;
  ofloat fblank = ObitMagicF();
  /* out = in1 (outer product) conjg(in2) */
  if ((in1[0]!=fblank) && (in2[0]!=fblank)) {
    out[0]  =  in1[0] * in2[0] + in1[1] * in2[1];  out[1]  =  in1[1] * in2[0] - in1[0] * in2[1];
    out[2]  =  in1[0] * in2[2] + in1[1] * in2[3];  out[3]  =  in1[1] * in2[2] - in1[0] * in2[3];
    out[4]  =  in1[2] * in2[0] + in1[3] * in2[1];  out[5]  =  in1[3] * in2[0] - in1[2] * in2[1];
    out[6]  =  in1[2] * in2[2] + in1[3] * in2[3];  out[7]  =  in1[3] * in2[2] - in1[2] * in2[3];
    
    out[8]  =  in1[0] * in2[4] + in1[1] * in2[5];  out[9]  =  in1[1] * in2[4] - in1[0] * in2[5];
    out[10] =  in1[0] * in2[6] + in1[1] * in2[7];  out[11] =  in1[1] * in2[6] - in1[0] * in2[7];
    out[12] =  in1[2] * in2[4] + in1[3] * in2[5];  out[13] =  in1[3] * in2[4] - in1[2] * in2[5];
    out[14] =  in1[2] * in2[6] + in1[3] * in2[7];  out[15] =  in1[3] * in2[6] - in1[2] * in2[7];
    
    out[16] =  in1[4] * in2[0] + in1[5] * in2[1];  out[17] =  in1[5] * in2[0] - in1[4] * in2[1];
    out[18] =  in1[4] * in2[2] + in1[5] * in2[3];  out[19] =  in1[5] * in2[2] - in1[4] * in2[3];
    out[20] =  in1[6] * in2[0] + in1[7] * in2[1];  out[21] =  in1[7] * in2[0] - in1[6] * in2[1];
    out[22] =  in1[6] * in2[2] + in1[7] * in2[3];  out[23] =  in1[7] * in2[2] - in1[6] * in2[3];
    
    out[24] =  in1[4] * in2[4] + in1[5] * in2[5];  out[25] =  in1[5] * in2[4] - in1[4] * in2[5];
    out[26] =  in1[4] * in2[6] + in1[5] * in2[7];  out[27] =  in1[5] * in2[6] - in1[4] * in2[7];
    out[28] =  in1[6] * in2[4] + in1[7] * in2[5];  out[29] =  in1[7] * in2[4] - in1[6] * in2[5];
    out[30] =  in1[6] * in2[6] + in1[7] * in2[7];  out[31] =  in1[7] * in2[6] - in1[6] * in2[7];
  } else {  /* One bad - blank */
    for (i=0; i<32; i++) out[i] = fblank;
  }
  } /* end  MatxOuter */

