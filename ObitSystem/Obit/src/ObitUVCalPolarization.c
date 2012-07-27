/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2012                                          */
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

#include "ObitUVCalPolarizationDef.h"
#include "ObitUVCalPolarization.h"
#include "ObitTablePD.h"
#include "ObitSinCos.h"
#include "ObitComplex.h"
#include "gsl/gsl_linalg.h"

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
  me->numIF     = desc->inaxes[desc->jlocif];
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
  olong itemp, SubA, SourID, FreqID, iChan, limit, index, jndex, loff;
  olong i, iif, ifreq, nch, ipol, ioff, joff, koff, jrl, jlr, ia1, ia2, ifoff;
  olong j, voff[4], choff, chdelta;
  gboolean wflag, someOK;
  ofloat xtemp[8], ytemp[8], dtemp[8], Lambda2, fblank = ObitMagicF();
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

  /* DEBUG unit poln vis */
  for (j=0; j<8; j++) xtemp[j] = 0.0;  xtemp[0] =  xtemp[6] = 1.0;

  /* Subarray number in data */
  itemp = (olong)RP[desc->ilocb];
  SubA = 1 + (olong)(100.0*(RP[desc->ilocb] -(ofloat)itemp) + 0.1);
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
      dtemp[2] = visIn[joff+6]; dtemp[3] = visIn[joff+7]; 
      dtemp[4] = visIn[joff+9]; dtemp[5] = visIn[joff+10]; 
      dtemp[6] = visIn[joff+3]; dtemp[7] = visIn[joff+4]; 
      
      /* Now apply calibration by multiplying the inverse Mueller matrix by 
	 the data vector. */
      j = 0;
      if (me->perChan) jndex = (ifoff + choff);
      else             jndex = ifoff;
      limit = 8 * in->numStok;
      MatxVec4Mult(&me->PolCal[jndex], dtemp, ytemp);
      /*MatxVec4Mult(&me->PolCal[jndex], xtemp, dtemp); DEBUG */
      
      index = 0;
      j = 0;
      limit = in->numStok;
      /* Reorder */
      visIn[joff+0] = ytemp[0];  visIn[joff+1]  = ytemp[1]; 
      visIn[joff+3] = ytemp[6];  visIn[joff+4]  = ytemp[7]; 
      visIn[joff+6] = ytemp[2];  visIn[joff+7]  = ytemp[3]; 
      visIn[joff+9] = ytemp[4];  visIn[joff+10] = ytemp[5]; 
      
      /* Done if in->numStok < 4) */
      if (in->numStok < 4) continue;
      
      /* parallactic angle  - not for 'ori-eli' and 'vlbi' */
      if ((me->polType==OBIT_UVPoln_ELORI) || (me->polType==OBIT_UVPoln_VLBI)) {
	gr = 1.0;
	gi = 0.0;
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
 * Multiplies by parallactic angle correction
 * Elipticity/Orientation model, one matrix per IF or channel/IF.
 * Input uses values in Jones arrays which are not a Jones matrix itself
 * Matrix assumes data order RR,RL,LR,LL
 * \param in      Polarization Object.
 * \param iant1   First antenna of baseline
 * \param iant2   Second antenna of baseline
 * \param iChan   Central channel to use
 * \param err     Error Stack.
 */
static void OriPol(ObitUVCalPolarizationS *in, olong iant1, olong iant2, olong iChan, ObitErr *err)
{
  olong iif, ia1, ia2, ich, bChan, eChan, kndx, jndx, nch;
  ofloat *M;
  ObitPolCalList *PCal = in->PCal;
  int ret, signum, ii, jj, kk, n;
  gsl_complex val;
  gsl_matrix_complex *inMatrix=NULL, *outMatrix=NULL;
  gsl_permutation *p=NULL;
  ofloat  PD;

  ocomplex PRref, PLref, PPRL, PPLR, PA1, PA2, PA1c, PA2c, ct1;
  ocomplex VRR, VRL, VLR, VLL;
  ocomplex RS1, RD1, LS1, LD1, RS2, RD2, LS2, LD2;
  ofloat root2, chi1, chi2;

  /* GSL work matrices */
  n = 4;
  p = gsl_permutation_alloc(n);
  inMatrix  = gsl_matrix_complex_alloc(n, n);
  outMatrix = gsl_matrix_complex_alloc(n, n);

  /* Per Channels? */
  if (in->perChan) {bChan = in->bChan; eChan = in->eChan; }
  else             {bChan = 1;         eChan = bChan; }


  /* WHERE IS THIS FROM??? */
  PD = 0.0;

  /* Parallactic angle correction */
  ia1 = iant1 - 1;
  ia2 = iant2 - 1;
  chi1 = in->curPA[ia1];
  chi2 = in->curPA[ia2];
  COMPLEX_EXP (PA1, -2*chi1);
  COMPLEX_EXP (PA2, -2*chi2);
  COMPLEX_CONJUGATE (PA1c, PA1);
  COMPLEX_CONJUGATE (PA2c, PA2);
  root2 = 1.0 / sqrt(2.0);
  /* loop thru IFs */
  jndx = 0; /* Offset in PolCal to beginning of IF matrix */
  M    = in->PolCal;
  for (iif=  in->bIF; iif<=in->eIF; iif++) { /* loop 400 */
 
    /* Loop over channel */
    nch  = in->numChan;
    for (ich=bChan; ich<=eChan; ich++) {
      kndx = 4*((iif-1)*nch + (ich-1));

      /* First antenna paramaters - values saved */
      COMPLEX_SET (RS1, in->Jones[ia1][jndx+0], in->Jones[ia1][jndx+1]);
      COMPLEX_SET (RD1, in->Jones[ia1][jndx+2], in->Jones[ia1][jndx+3]);
      COMPLEX_SET (LS1, in->Jones[ia1][jndx+4], in->Jones[ia1][jndx+5]);
      COMPLEX_SET (LD1, in->Jones[ia1][jndx+6], in->Jones[ia1][jndx+7]);

      /* conjugate of second antenna - values saved */
      COMPLEX_SET (RS2, in->Jones[ia2][jndx+0], -in->Jones[ia2][jndx+1]);
      COMPLEX_SET (RD2, in->Jones[ia2][jndx+2], -in->Jones[ia2][jndx+3]);
      COMPLEX_SET (LS2, in->Jones[ia2][jndx+4], -in->Jones[ia2][jndx+5]);
      COMPLEX_SET (LD2, in->Jones[ia2][jndx+6], -in->Jones[ia2][jndx+7]);

      /* Reference antenna phase terms */
      COMPLEX_EXP (PRref,  PCal->ANlist[PCal->polRefAnt-1][kndx+1]);
      COMPLEX_EXP (PLref, -PCal->ANlist[PCal->polRefAnt-1][kndx+3]+PD);
      COMPLEX_CONJUGATE (ct1, PLref);
      COMPLEX_MUL2 (PPRL, PRref, ct1);
      COMPLEX_CONJUGATE (ct1, PRref);
      COMPLEX_MUL2 (PPLR, PLref, ct1);

      /* Terms of Mueller matrix */
      /* VRR = RS1 * RS2        
  	       RS1 * RD2 * PA2c 
	       RD1 * RS2 * PA1  
	       RD1 * RD2 * PA1  * PA2c */
      COMPLEX_MUL2 (VRR, RS1,  RS2);
      GSL_SET_COMPLEX(&val, VRR.real, VRR.imag);
      gsl_matrix_complex_set (inMatrix, 0, 0, val);
      COMPLEX_MUL3 (VRR, RS1, RD2, PA2c);
      GSL_SET_COMPLEX(&val, VRR.real, VRR.imag);
      gsl_matrix_complex_set (inMatrix, 0, 1, val);
      COMPLEX_MUL3 (VRR, RD1,  RS2,  PA1);
      GSL_SET_COMPLEX(&val, VRR.real, VRR.imag);
      gsl_matrix_complex_set (inMatrix, 0, 2, val);
      COMPLEX_MUL4 (VRR, RD1,  RD2, PA1,  PA2c);
      GSL_SET_COMPLEX(&val, VRR.real, VRR.imag);
      gsl_matrix_complex_set (inMatrix, 0, 3, val);
      
      /* VLL = LS1 * LS2 * PA1c * PA2	
	       LS1 * LD2 * PA1c
	       LD1 * LS2 * PA2 
	       LD1 * LD2 */
      COMPLEX_MUL4 (VLL, LS1,  LS2, PA1c,  PA2);
      GSL_SET_COMPLEX(&val, VLL.real, VLL.imag);
      gsl_matrix_complex_set (inMatrix, 3, 0, val);
      COMPLEX_MUL3 (VLL, LS1,  LD2,  PA1c);
      GSL_SET_COMPLEX(&val, VLL.real, VLL.imag);
      gsl_matrix_complex_set (inMatrix, 3, 1, val);
      COMPLEX_MUL3 (VLL, LD1, LS2,  PA2);
      GSL_SET_COMPLEX(&val, VLL.real, VLL.imag);
      gsl_matrix_complex_set (inMatrix, 3, 2, val);
      COMPLEX_MUL2 (VLL, LD1,  LD2);
      GSL_SET_COMPLEX(&val, VLL.real, VLL.imag);
      gsl_matrix_complex_set (inMatrix, 3, 3, val);
      
      /* VRL = PPRL * RS1 * LS2 * PA2
   	       PPRL * RS1 * LD2 
	       PPRL * RD1 * LS2 * PA1 * PA2
	       PPRL * RD1 * LD2 * PA1 */
      COMPLEX_MUL2 (ct1, PPRL, RS1);
      COMPLEX_MUL3 (VRL, ct1,  LS2,  PA2);
      GSL_SET_COMPLEX(&val, VRL.real, VRL.imag);
      gsl_matrix_complex_set (inMatrix, 1, 0, val);
      COMPLEX_MUL2 (VRL, ct1,  LD2);
      GSL_SET_COMPLEX(&val, VRL.real, VRL.imag);
      gsl_matrix_complex_set (inMatrix, 1, 1, val);
      COMPLEX_MUL3 (ct1, PPRL, RD1, PA1);
      COMPLEX_MUL3 (VRL, ct1, LS2,  PA2);
      GSL_SET_COMPLEX(&val, VRL.real, VRL.imag);
      gsl_matrix_complex_set (inMatrix, 1, 2, val);
      COMPLEX_MUL2 (VRL, ct1, LD2);
      GSL_SET_COMPLEX(&val, VRL.real, VRL.imag);
      gsl_matrix_complex_set (inMatrix, 1, 3, val);
      
      /* VLR = PPLR * LS1 * RS2 * PA1c
	       PPLR * LS1 * RD2 * PA1c * PA2c
	       PPLR * LD1 * RS2
	       PPLR * LD1 * RD2 * PA2c */
      COMPLEX_MUL3 (ct1, PPLR, LS1,  PA1c);
      COMPLEX_MUL2 (VLR, ct1,  RS2);
      GSL_SET_COMPLEX(&val, VLR.real, VLR.imag);
      gsl_matrix_complex_set (inMatrix, 2, 0, val);
      COMPLEX_MUL3 (VLR, ct1, RD2, PA2c);
      GSL_SET_COMPLEX(&val, VLR.real, VLR.imag);
      gsl_matrix_complex_set (inMatrix, 2, 1, val);
      COMPLEX_MUL2 (ct1, PPLR, LD1);
      COMPLEX_MUL2 (VLR, ct1,  RS2);
      GSL_SET_COMPLEX(&val, VLR.real, VLR.imag);
      gsl_matrix_complex_set (inMatrix, 2, 2, val);
      COMPLEX_MUL3 (VLR,  ct1, RD2, PA2c);
      GSL_SET_COMPLEX(&val, VLR.real, VLR.imag);
      gsl_matrix_complex_set (inMatrix, 2, 3, val);
      
      /* Invert Mueller matrix  */
      ret = gsl_linalg_complex_LU_decomp (inMatrix, p, &signum);
      ret = gsl_linalg_complex_LU_invert (inMatrix, p, outMatrix);

      /* get values */
      kk = 0;
      for (ii=0; ii<4; ii++) {
	for (jj=0; jj<4; jj++) {
	  val = gsl_matrix_complex_get (outMatrix, ii, jj);
	  /*val = gsl_matrix_complex_get (inMatrix, ii, jj);*/
	  M[kk++]  = GSL_REAL(val);
	  M[kk++]  = GSL_IMAG(val);
	}
      }

      M += 32;    /* Offset in PolCal to beginning of channel/IF matrix */
      jndx +=  8; /* Offset in Jones to beginning of  channel/IF  */
    } /* end channel loop */
  } /* end loop over IF */
  if (inMatrix)  gsl_matrix_complex_free(inMatrix);
  if (outMatrix) gsl_matrix_complex_free(outMatrix);
  if (p)         gsl_permutation_free(p);
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
      
      /* iant1 * conjg(iant2), for order RR,RL,LR,LL */
      in->PolCal[jndex]    =  in->Jones[ia1][ioff+0] * in->Jones[ia2][ioff+0] + in->Jones[ia1][ioff+1] * in->Jones[ia2][ioff+1];
      in->PolCal[jndex+1]  =  in->Jones[ia1][ioff+0] * in->Jones[ia2][ioff+1] - in->Jones[ia1][ioff+1] * in->Jones[ia2][ioff+0];
      in->PolCal[jndex+2]  =  in->Jones[ia1][ioff+0] * in->Jones[ia2][ioff+2] + in->Jones[ia1][ioff+1] * in->Jones[ia2][ioff+3];
      in->PolCal[jndex+3]  =  in->Jones[ia1][ioff+0] * in->Jones[ia2][ioff+3] - in->Jones[ia1][ioff+1] * in->Jones[ia2][ioff+2];
      in->PolCal[jndex+4]  =  in->Jones[ia1][ioff+2] * in->Jones[ia2][ioff+0] + in->Jones[ia1][ioff+3] * in->Jones[ia2][ioff+1];
      in->PolCal[jndex+5]  =  in->Jones[ia1][ioff+2] * in->Jones[ia2][ioff+1] - in->Jones[ia1][ioff+3] * in->Jones[ia2][ioff+0];
      in->PolCal[jndex+6]  =  in->Jones[ia1][ioff+2] * in->Jones[ia2][ioff+2] + in->Jones[ia1][ioff+3] * in->Jones[ia2][ioff+3];
      in->PolCal[jndex+7]  =  in->Jones[ia1][ioff+2] * in->Jones[ia2][ioff+3] - in->Jones[ia1][ioff+3] * in->Jones[ia2][ioff+2];
      
      in->PolCal[jndex+8]  =  in->Jones[ia1][ioff+0] * in->Jones[ia2][ioff+4] + in->Jones[ia1][ioff+1] * in->Jones[ia2][ioff+5];
      in->PolCal[jndex+9]  =  in->Jones[ia1][ioff+0] * in->Jones[ia2][ioff+5] - in->Jones[ia1][ioff+1] * in->Jones[ia2][ioff+4];
      in->PolCal[jndex+10] =  in->Jones[ia1][ioff+0] * in->Jones[ia2][ioff+6] + in->Jones[ia1][ioff+1] * in->Jones[ia2][ioff+7];
      in->PolCal[jndex+11] =  in->Jones[ia1][ioff+0] * in->Jones[ia2][ioff+7] - in->Jones[ia1][ioff+1] * in->Jones[ia2][ioff+6];
      in->PolCal[jndex+12] =  in->Jones[ia1][ioff+2] * in->Jones[ia2][ioff+4] + in->Jones[ia1][ioff+3] * in->Jones[ia2][ioff+5];
      in->PolCal[jndex+13] =  in->Jones[ia1][ioff+2] * in->Jones[ia2][ioff+5] - in->Jones[ia1][ioff+3] * in->Jones[ia2][ioff+4];
      in->PolCal[jndex+14] =  in->Jones[ia1][ioff+2] * in->Jones[ia2][ioff+6] + in->Jones[ia1][ioff+3] * in->Jones[ia2][ioff+7];
      in->PolCal[jndex+15] =  in->Jones[ia1][ioff+2] * in->Jones[ia2][ioff+7] - in->Jones[ia1][ioff+3] * in->Jones[ia2][ioff+6];
      
      in->PolCal[jndex+16] =  in->Jones[ia1][ioff+4] * in->Jones[ia2][ioff+0] + in->Jones[ia1][ioff+5] * in->Jones[ia2][ioff+1];
      in->PolCal[jndex+17] =  in->Jones[ia1][ioff+4] * in->Jones[ia2][ioff+1] - in->Jones[ia1][ioff+5] * in->Jones[ia2][ioff+0];
      in->PolCal[jndex+18] =  in->Jones[ia1][ioff+4] * in->Jones[ia2][ioff+2] + in->Jones[ia1][ioff+5] * in->Jones[ia2][ioff+3];
      in->PolCal[jndex+19] =  in->Jones[ia1][ioff+4] * in->Jones[ia2][ioff+3] - in->Jones[ia1][ioff+5] * in->Jones[ia2][ioff+2];
      in->PolCal[jndex+20] =  in->Jones[ia1][ioff+6] * in->Jones[ia2][ioff+0] + in->Jones[ia1][ioff+7] * in->Jones[ia2][ioff+1];
      in->PolCal[jndex+21] =  in->Jones[ia1][ioff+6] * in->Jones[ia2][ioff+1] - in->Jones[ia1][ioff+7] * in->Jones[ia2][ioff+1];
      in->PolCal[jndex+22] =  in->Jones[ia1][ioff+6] * in->Jones[ia2][ioff+2] + in->Jones[ia1][ioff+7] * in->Jones[ia2][ioff+3];
      in->PolCal[jndex+23] =  in->Jones[ia1][ioff+6] * in->Jones[ia2][ioff+3] - in->Jones[ia1][ioff+7] * in->Jones[ia2][ioff+2];
      
      in->PolCal[jndex+24] =  in->Jones[ia1][ioff+4] * in->Jones[ia2][ioff+4] + in->Jones[ia1][ioff+5] * in->Jones[ia2][ioff+5];
      in->PolCal[jndex+25] =  in->Jones[ia1][ioff+4] * in->Jones[ia2][ioff+5] - in->Jones[ia1][ioff+5] * in->Jones[ia2][ioff+4];
      in->PolCal[jndex+26] =  in->Jones[ia1][ioff+4] * in->Jones[ia2][ioff+6] + in->Jones[ia1][ioff+5] * in->Jones[ia2][ioff+7];
      in->PolCal[jndex+27] =  in->Jones[ia1][ioff+4] * in->Jones[ia2][ioff+7] - in->Jones[ia1][ioff+5] * in->Jones[ia2][ioff+6];
      in->PolCal[jndex+28] =  in->Jones[ia1][ioff+6] * in->Jones[ia2][ioff+4] + in->Jones[ia1][ioff+7] * in->Jones[ia2][ioff+5];
      in->PolCal[jndex+29] =  in->Jones[ia1][ioff+6] * in->Jones[ia2][ioff+5] - in->Jones[ia1][ioff+7] * in->Jones[ia2][ioff+4];
      in->PolCal[jndex+30] =  in->Jones[ia1][ioff+6] * in->Jones[ia2][ioff+6] + in->Jones[ia1][ioff+7] * in->Jones[ia2][ioff+7];
      in->PolCal[jndex+31] =  in->Jones[ia1][ioff+6] * in->Jones[ia2][ioff+7] - in->Jones[ia1][ioff+7] * in->Jones[ia2][ioff+6];
      
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
  ofloat Jones[8], Dr[2]={0.0,0.0}, Dl[2]={0.0,0.0}, Det[2], d;
  ofloat elp_r, elp_l, ori_r, ori_l, angle[4], sina[4], cosa[4];
  ofloat root2, fblank = ObitMagicF();
  ofloat rotate=0.0, crot, srot, temp[8];
  olong iif, jndx, loff, refAnt, i, SubA;
  gboolean doJones;

  SubA = MAX (1, in->curSubA);
  refAnt = Ant->polRefAnt;
  root2  = 1.0 / sqrt(2.0);

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
      ori_r = Ant->ANlist[refAnt-1]->FeedAPCal[2*iif+1];
      elp_l = Ant->ANlist[iant-1]->FeedBPCal[2*iif+0];
      ori_l = Ant->ANlist[refAnt-1]->FeedBPCal[2*iif+1] + Ant->RLPhaseDiff[iif];
      
      /* Sines/cosines */
      angle[0] = elp_r; angle[1] = 2*ori_r; angle[2] = elp_l; angle[3] = -2*ori_l;
      ObitSinCosVec(4, angle, sina, cosa);
      
      /* Jones matrix: 
	 (root2*(cos(Er)+sin(Er))+i0)              (root2*(cos(Er)-sin(Er)), i0)*exp(i2*Or)
	 (root2*(cos(El)+sin(El))+i0)*exp(-i2*Ol)  (root2*(cos(El)-sin(El)), i0)
      */
      in->Jones[iant-1][jndx+0] = root2*(cosa[0] + sina[0]);
      in->Jones[iant-1][jndx+1] = 0.0;
      in->Jones[iant-1][jndx+2] = root2*(cosa[0] - sina[0]) * cosa[1];
      in->Jones[iant-1][jndx+3] = root2*(cosa[0] - sina[0]) * sina[1];
      in->Jones[iant-1][jndx+4] = root2*(cosa[2] + sina[2]) * cosa[3];
      in->Jones[iant-1][jndx+5] = root2*(cosa[2] + sina[2]) * sina[3];
      in->Jones[iant-1][jndx+6] = root2*(cosa[2] - sina[2]);
      in->Jones[iant-1][jndx+7] = 0.0;

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
      
      /* invert matrix */
      in->Jones[iant-1][jndx+6] =   Jones[0] * Det[0] - Jones[1] * Det[1];
      in->Jones[iant-1][jndx+7] =   Jones[0] * Det[1] + Jones[1] * Det[0];
      in->Jones[iant-1][jndx+2] = -(Jones[2] * Det[0] - Jones[3] * Det[1]);
      in->Jones[iant-1][jndx+3] = -(Jones[2] * Det[1] + Jones[3] * Det[0]);
      in->Jones[iant-1][jndx+4] = -(Jones[4] * Det[0] - Jones[5] * Det[1]);
      in->Jones[iant-1][jndx+5] = -(Jones[4] * Det[1] + Jones[5] * Det[0]);
      in->Jones[iant-1][jndx+0] =   Jones[6] * Det[0] - Jones[7] * Det[1];
      in->Jones[iant-1][jndx+1] =   Jones[6] * Det[1] + Jones[7] * Det[0];
    } /* end doJones */
    
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
  ofloat elp_r, elp_l, ori_r, ori_l, angle[4], sina[4], cosa[4];
  ofloat root2, fblank = ObitMagicF();
  ofloat rotate=0.0, crot, srot, temp[8];
  olong iif, jndx, loff, refAnt, i, SubA, ich, nch, kndx, ia;
  gboolean doJones;
  ObitPolCalList *PCal = in->PCal;
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
  nch    = in->numChan;
  ia     = iant - 1;
  root2  = 1.0 / sqrt(2.0);

  /* Loop over IFs (index 0 rel) */
  jndx = 0;
  for (iif=in->bIF-1; iif<=in->eIF-1; iif++) {

    /* Loop over channel */
    for (ich=0; ich<nch; ich++) {
      kndx = 4*(iif*nch + ich + (in->bChan-1));
      
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

	/* Sines/cosines */
	angle[0] = elp_r; angle[1] = 2*ori_r; angle[2] = elp_l; angle[3] = -2*ori_l;
	ObitSinCosVec(4, angle, sina, cosa);

	/* Jones matrix: 
	   (root2*(cos(Er)+sin(Er))+i0)              (root2*(cos(Er)-sin(Er)), i0)*exp(i2*Or)
	   (root2*(cos(El)+sin(El))+i0)*exp(-i2*Ol)  (root2*(cos(El)-sin(El)), i0)
	 */
 	in->Jones[ia][jndx+0] = root2*(cosa[0] + sina[0]);
	in->Jones[ia][jndx+1] = 0.0;
	in->Jones[ia][jndx+2] = root2*(cosa[0] - sina[0]) * cosa[1];
	in->Jones[ia][jndx+3] = root2*(cosa[0] - sina[0]) * sina[1];
	in->Jones[ia][jndx+4] = root2*(cosa[2] + sina[2]) * cosa[3];
 	in->Jones[ia][jndx+5] = root2*(cosa[2] + sina[2]) * sina[3];
	in->Jones[ia][jndx+6] = root2*(cosa[2] - sina[2]);
	in->Jones[ia][jndx+7] = 0.0;

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
      
      /* Need to compute Jones matric ain inverse determinant? */
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
	
	/* invert matrix */
	in->Jones[ia][jndx+6] =   Jones[0] * Det[0] - Jones[1] * Det[1];
	in->Jones[ia][jndx+7] =   Jones[0] * Det[1] + Jones[1] * Det[0];
	in->Jones[ia][jndx+2] = -(Jones[2] * Det[0] - Jones[3] * Det[1]);
	in->Jones[ia][jndx+3] = -(Jones[2] * Det[1] + Jones[3] * Det[0]);
	in->Jones[ia][jndx+4] = -(Jones[4] * Det[0] - Jones[5] * Det[1]);
	in->Jones[ia][jndx+5] = -(Jones[4] * Det[1] + Jones[5] * Det[0]);
	in->Jones[ia][jndx+0] =   Jones[6] * Det[0] - Jones[7] * Det[1];
	in->Jones[ia][jndx+1] =   Jones[6] * Det[1] + Jones[7] * Det[0];
     } /* end doJones */
     
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

