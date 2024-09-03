/* $Id$ */
/* To do 
- implement Ionospheric Faraday rotation calibration
  o ObitUVCalJones - cal for linear feeds, update table
 */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2024                                               */
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

#include "ObitMatx.h"
#include "ObitComplex.h"
#include "ObitUVCalJonesDef.h"
#include "ObitUVCalJones.h"
#include "ObitTableJI.h"
#include "ObitTableJT.h"
#include "ObitTableUtil.h"
#include "ObitUVDesc.h"
#include "ObitSinCos.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVCalJones.c
 * ObitUVCal utilities for applying Jones matrix calibration to uv data.
 */

/*-------------------Private function prototypes-------------------*/
/** Private:  Create structure for Jones matrix calibration . */
static ObitUVCalJonesS* newObitUVCalJonesS (ObitUVCal *in);

/** Private: Update calibration arrays. */
static void ObitUVCalJonesUpdate (ObitUVCalJonesS *in, ObitUVCalCalibrateS *cal,
				  ObitUVSel *sel, ofloat time, olong SourID, olong FreqID,
				  ObitErr *err);

/** Private:  Read calibration for a new time into the internal arrays. */
static void ObitUVCalJonesNewTime (ObitUVCalJonesS *in, ObitUVCalCalibrateS *cal,
				   ObitUVSel *sel, ofloat time,  olong SourID, olong FreqID,
				   ObitErr *err);

/*----------------------Public functions---------------------------*/
/**
 * Initialize structures for Jones matrix Calibration .
 * \param in   Calibration Object.
 * \param sel  Data selector.
 * \param desc Data descriptor.
 * \param err  ObitError stack.
 */
void ObitUVCalJonesInit (ObitUVCal *in, ObitUVSel *sel, ObitUVDesc *desc, 
			 ObitErr *err)
{
  /*  ObitIOCode retCode;*/
  ObitUVCalJonesS *me;
  ofloat t1=0.0, t2, t3;
  olong i, size, nant, naxis[2] = {2,2};
  gchar *colName="TIME    ";
  gchar *routine = "ObitUVCalJonesInit";

  /* Is Jones calibration requested? */
  if (!sel->doJones) return;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVCalIsA(in));

  /* Make sure the Jones calibration Table exists */
  if (!in->JTable) {
    Obit_log_error(err, OBIT_Error, "%s: Jones Calibration Table not provided", routine);
    return;
  }

  in->doPol = TRUE;

  in->JonesCal = newObitUVCalJonesS(in);

  /* pointer to calibration structure */
  me = in->JonesCal;

  /* Copy Selector information */
  me->bChan       = sel->startChann;
  me->eChan       = sel->startChann + sel->numberChann - 1;
  me->bIF         = sel->startIF;
  me->eIF         = sel->startIF + sel->numberIF - 1;
  me->SubA        = sel->SubA;
  me->FreqID      = sel->FreqID;
  me->useJI       = FALSE;
  me->useJT       = FALSE;

  /* Copy Cal information */
  me->numSubA     = in->numSubA;

  /* Copy descriptor information */
  me->numAnt    = desc->maxAnt;
  me->numSubA   = desc->numSubA;
  if (desc->jlocif>=0) me->numIF = desc->inaxes[desc->jlocif];
  else                 me->numIF = 1;
  me->numChan   = desc->inaxes[desc->jlocf];

  me->curPA    = g_realloc(me->curPA,    me->numAnt*sizeof(ofloat));
  me->curCosPA = g_realloc(me->curCosPA, me->numAnt*sizeof(ofloat));
  me->curSinPA = g_realloc(me->curSinPA, me->numAnt*sizeof(ofloat));
  size = me->numAnt * (in->eIF-in->bIF+1)*(in->eChan-in->bChan+1);
  me->Jones    = g_realloc(me->Jones,    size*sizeof(ObitMatx*));
  for (i=0; i<size; i++) me->Jones[i]  = ObitMatxCreate(OBIT_Complex, 2, naxis);
  me->TempMatx  = ObitMatxCreate(OBIT_Complex, 2, naxis);
  me->TempMatx2 = ObitMatxCreate(OBIT_Complex, 2, naxis);

  /* Init time, source... */
  me->curTime       = -1.0e20;
  me->curSourID     = -1;
  me->curSubA       = -1;
  me->CalTime       = -1.0e5;
  me->PriorCalTime  = -1.0e5;
  me->FollowCalTime = -1.0e5;
  me->LastRowRead   = -10;
  me->numRow         = 0;

  /* Init Sine/cosine function */
  ObitSinCosCalc(t1, &t2, &t3);

  /* Are the data from circular or linear feeds */
  me->circFeed = desc->crval[desc->jlocs] > -4.0;

  /* Is this MeerKAT data? */
  me->isMeerKAT = !strncmp(in->antennaLists[0]->ArrName, "MeerKAT", 7);

  /* Cal Table AIPS JI for single source */
  if (desc->ilocsu<0) {  /* Source ID random parameter? */
    me->JITable = ObitTableRef(in->JTable);
    ObitTableUtilSort((ObitTable*)(me->JITable), colName, FALSE, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    nant = ((ObitTableJI*)me->JITable)->numAnt;
    ObitTableJIOpen ((ObitTableJI*)(me->JITable), OBIT_IO_ReadOnly, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    ((ObitTableJI*)me->JITable)->numAnt = 
      MAX(((ObitTableJI*)me->JITable)->numAnt, nant);  /* Don't let it get smaller */
    me->JITableRow = (Obit*)newObitTableJIRow((ObitTableJI*)(me->JITable));
    me->numRow = ((ObitTableJI*)me->JITable)->myDesc->nrow;
    me->useJI = TRUE;
    /* Check number of antennas */
    Obit_return_if_fail((me->numAnt<=((ObitTableJI*)me->JITable)->numAnt), err,
			"%s: JI Table different number of antennas < data %d, %d",
			routine, ((ObitTableJI*)me->JITable)->numAnt, me->numAnt);
    /* Tell if prtLv>=1 */
    if (err->prtLv>=1) {
      Obit_log_error(err, OBIT_InfoErr, "Calibrating with Jones Table JI %d",sel->JonesVersion);
      ObitErrLog(err);
    }
  } else {   /* AIPS JT for multisource */
    me->JTTable = ObitTableRef(in->JTable);
    ObitTableUtilSort((ObitTable*)(me->JTTable), colName, FALSE, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    nant = ((ObitTableJT*)me->JTTable)->numAnt;
    ObitTableJTOpen ((ObitTableJT*)(me->JTTable), OBIT_IO_ReadOnly, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    ((ObitTableJT*)me->JTTable)->numAnt = 
      MAX(((ObitTableJT*)me->JTTable)->numAnt, nant);  /* Don't let it get smaller */
    me->JTTableRow = (Obit*)newObitTableJTRow((ObitTableJT*)(me->JTTable));
    me->numRow = ((ObitTableJT*)me->JTTable)->myDesc->nrow;
    me->useJT = TRUE;
    /* Check number of antennas */
    Obit_return_if_fail((me->numAnt<=((ObitTableJT*)me->JTTable)->numAnt), err,
			"%s: JT Table different number of antennas < data %d, %d",
			routine, ((ObitTableJT*)me->JTTable)->numAnt, me->numAnt);
    /* Tell if prtLv>=1 */
    if (err->prtLv>=1) {
      Obit_log_error(err, OBIT_InfoErr, "Calibrating with Jones Table JT %d",sel->JonesVersion);
      ObitErrLog(err);
    }
 } /* end AIPS JT table */

} /*  end ObitUVCalJonesInit */

/**
 * Calibrate data using Jones matrices.
 *
 * \param in    Jones Matrix Calibration Object.
 * \param time  Time of datum
 * \param ant1  first antenna number of baseline
 * \param ant2  second antanna of baseline.
 * \param RP    Random parameters array.
 * \param visIn 1 visibility as an array of floats
 * \param err   ObitError stack.
 */
void ObitUVCalJones (ObitUVCal *in, float time, olong ant1, olong ant2, 
		     ofloat *RP, ofloat *visIn, ObitErr *err)
{
  olong SubA, loff, i, iif, ifreq, ioff, joff=0, jrl, jlr, ia1, ia2;
  olong  FreqID, SourID, it1, it2;
  /* olong FreqID, SourID, ifoff, choff;*/
  olong indx1, indx2, jndx1, jndx2, kndx1, kndx2, nChIF;
  /*olong  voff[4], */
  gboolean wflag, someOK, flagged=FALSE, doIFR;
  ofloat Lambda2, gr, gi, gr1, gi1, tr, ti;
  ObitUVCalJonesS *me;
  ObitUVCalCalibrateS *cal;
  ObitUVDesc *desc;
  gchar *routine="ObitUVCalJones";

  if (err->error) return;

  /* local pointers for structures */
  me   = in->JonesCal;
  desc = in->myDesc;
  cal  = in->ampPhaseCal;

  /* Make sure some good data */
  someOK = FALSE;
  for (i=0; i<desc->ncorr; i++) someOK = someOK || (visIn[i*3+2]>0.0);
  if (!someOK) return;
  
  /* Data Freq id */
  if (desc->ilocfq >= 0) FreqID = RP[desc->ilocfq] + 0.1;
  else  FreqID = 0;
  
  /* Source ID */
  if (desc->ilocsu >= 0) SourID = RP[desc->ilocsu] + 0.1;
  else SourID = 0;
  
  /* see if new time - update cal. */
  if (time > me->CalTime) {
    ObitUVCalJonesUpdate (me, cal, in->mySel, time, SourID, FreqID,  err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }
  
  /* Subarray number in data */
  ObitUVDescGetAnts(desc, RP, &it1, &it2, &SubA);
  SubA = MIN (SubA, in->numANTable);
  ia1  = ant1 - 1;
  ia2  = ant2 - 1;
  nChIF = (in->eIF-in->bIF+1) * (in->eChan-in->bChan+1); /* Number of channels*IFs */

  /* loop thru IF */
  kndx1 = ia1 * nChIF; kndx2 = ia2 * nChIF; 
  for (iif=me->bIF; iif<=me->eIF; iif++) { /* loop thru IFs */
    ioff = (iif-1) * desc->incif;
    jndx1 = kndx1+(iif-me->bIF)*(in->eChan-in->bChan+1);
    jndx2 = kndx2+(iif-me->bIF)*(in->eChan-in->bChan+1);
    
    /* loop thru channels */
    for (ifreq=me->bChan; ifreq<=me->eChan; ifreq++) {  /* channel loop */
      joff = ((ifreq-1) * desc->incf + ioff); /* Offset of RR (or XX) */
      indx1 = jndx1 + (ifreq-me->bChan); indx2 = jndx2 + (ifreq-me->bChan); 
      
      /* Any data flagged*/
      wflag = (visIn[joff+2]<=0.) || (visIn[joff+5]<=0.) || 
	(visIn[joff+8]<=0.) || (visIn[joff+11]<=0.);
      
      /* Either solution flagged? */
      flagged = ObitMatxIsBlank(me->Jones[indx1]) || ObitMatxIsBlank(me->Jones[indx2]);
      flagged = flagged || wflag;
      
      if (!flagged) {
	/* Save in work Matrix with reordering */
	ObitMatxSet2C(me->TempMatx, visIn[joff+0], visIn[joff+1], visIn[joff+6], visIn[joff+7],
		      visIn[joff+9], visIn[joff+10], visIn[joff+3], visIn[joff+4]);
	/* Problem with XY/YX swapped??? 
	ObitMatxSet2C(me->TempMatx, visIn[joff+0], visIn[joff+1], visIn[joff+9], visIn[joff+10],
		      visIn[joff+6], visIn[joff+7], visIn[joff+3], visIn[joff+4]);*/
	
	/* Calibrate */
	ObitMatxMult(me->Jones[indx1], me->TempMatx, me->TempMatx2);
	ObitMatxMultCT(me->TempMatx2, me->Jones[indx2], me->TempMatx);
	
	/* Extract */
	ObitMatxGet2C(me->TempMatx, &visIn[joff+0], &visIn[joff+1], &visIn[joff+6], &visIn[joff+7],
	  &visIn[joff+9], &visIn[joff+10], &visIn[joff+3], &visIn[joff+4]);
	/* Problem with XY/YX swapped???
 	ObitMatxGet2C(me->TempMatx, &visIn[joff+0], &visIn[joff+1], &visIn[joff+9], &visIn[joff+10],
	&visIn[joff+6], &visIn[joff+7], &visIn[joff+3], &visIn[joff+4]); */
      } else { /* end calibrate if not flagged */
	for (i=0; i<12; i++) visIn[i] = 0.0;   /* zero flagged data */
      }
    } /* end channel loop */
    
    
    /* IFR stuff in cal an ObitUVCalCalibrateS - THIS OUGHT TO BE DONE PER CHANNEL */
    /* REALLY NYI */ 
    doIFR = (cal!=NULL) && (cal->Lambda!=NULL) && (cal->Lambda[0]>0.0); /* Ionospheric cal? */
    /* correct RL,LR ionospheric Faraday rotation: */
    if (doIFR) {
      gr = 1.0; gi = 0.0;  /* Probably not needed */
      loff = (iif - 1) * cal->numLambda + ifreq - 1;
      Lambda2 = cal->Lambda[loff]*cal->Lambda[loff];
      gr1 = gr * cos (Lambda2 * (cal->IFR[ia1] + cal->IFR[ia2])) -  
	gi * sin (Lambda2 * (cal->IFR[ia1] + cal->IFR[ia2]));
      gi1 = gi * cos (Lambda2 * (cal->IFR[ia1] + cal->IFR[ia2])) +  
	gr * sin (Lambda2 * (cal->IFR[ia1] + cal->IFR[ia2]));
      
      /* Circular feeds - NEED ALSO FOR LINEAR*/
      if (me->circFeed) {
	/* Correct RL  */
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
      } else {  /* Linear feeds */
	g_error ("IFR not implemented for linear feeds - write me!");
      } /* end by feed type */
    } /* end loop over channels */
    jndx1++; jndx2++;  /* Increment IF indices */
  } /* end loop  over IF  */
} /* end ObitUVCalJones */

/**
 * Shutdown Jones Matrix calibration
 * Destroy structures.
 * \param in   Jones Matrix Object.
 * \param err  ObitError stack.
 */
void ObitUVCalJonesShutdown (ObitUVCal *in, ObitErr *err)
{
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVCalIsA(in));

  /* delete structure */
  in->JonesCal = ObitUVCalJonesSUnref(in->JonesCal);
} /*  end ObitUVCalJonesShutdown */

/**
 * Destroy structure for Jones matrix calibration.
 * \param in   Jones Matrix Object.
 * \return NULL
 */
ObitUVCalJonesS*
ObitUVCalJonesSUnref (ObitUVCalJonesS *in)
{
  olong i, num;

  if (in==NULL) return in;
  
  if (in->curPA)         {g_free(in->curPA);}         in->curPA = NULL;
  if (in->curCosPA)      {g_free(in->curCosPA);}      in->curCosPA = NULL;
  if (in->curSinPA)      {g_free(in->curSinPA);}      in->curSinPA = NULL;
  if (in->CalPrior)      {g_free(in->CalPrior);}      in->CalPrior = NULL;
  if (in->CalFollow)     {g_free(in->CalFollow);}     in->CalFollow = NULL;
  if (in->PriorAntTime)  {g_free(in->PriorAntTime);}  in->PriorAntTime = NULL;
  if (in->FollowAntTime) {g_free(in->FollowAntTime);} in->FollowAntTime = NULL;
  if (in->Jones) {
    num = in->numAnt * (in->eIF-in->bIF+1) * (in->eChan-in->bChan+1);
    for (i=0; i<num; i++) if (in->Jones[i]) ObitMatxUnref(in->Jones[i]);
    g_free(in->Jones); in->Jones = NULL;
  }
  /* Other Cal structures */
  in->TempMatx   = ObitMatxUnref(in->TempMatx);  in->TempMatx  = NULL;
  in->TempMatx2  = ObitMatxUnref(in->TempMatx2); in->TempMatx2 = NULL;
  if (in->useJI) {
    in->JITable    = ObitTableJIUnref((ObitTableJI*)in->JITable);
    in->JITableRow = ObitTableJIRowUnref((ObitTableJIRow*)in->JITableRow);
  }
  if (in->useJT) {
    in->JTTable    = ObitTableJTUnref((ObitTableJT*)in->JTTable);
    in->JTTableRow = ObitTableJTRowUnref((ObitTableJTRow*)in->JTTableRow);
  }
   
  /* basic structure */
   g_free (in);

  return NULL;
} /*  end ObitUVCalJonesSUnref */

/*---------------Private functions---------------------------*/

/**
 * Create structure for Jones Matrix Calibrate .
 * \param in   Jones Matrix Object.
 * \return newly created object.
 */
static ObitUVCalJonesS* newObitUVCalJonesS (ObitUVCal *in)
{
  ObitUVCalJonesS* out;

  out = g_malloc0(sizeof(ObitUVCalJonesS));
  out->useJI  = FALSE;
  out->useJT  = FALSE;
  out->numRow = 0;

  /* Init pointers */
  out->curPA        = NULL;
  out->curCosPA     = NULL;
  out->curSinPA     = NULL;
  out->Jones        = NULL;
  out->TempMatx     = NULL;
  out->TempMatx2    = NULL;
  out->JITable      = NULL;
  out->JITableRow   = NULL;
  out->JTTable      = NULL;
  out->JTTableRow   = NULL;
  out->PriorAntTime = NULL;
  out->FollowAntTime= NULL;
  out->CalPrior     = NULL;
  out->CalFollow    = NULL;
  return out;
} /*  end newObitUVCalJonesS */

/**
 * Update ObitUVCalJones calibration tables for time time.
 * The current JI/JT table is interpolated between the previous and following
 * sets of solutions.
 * If a new set of entries is needed from the JI/JT table they are read.
 * A calibration array entry (prior or following) consists of:
 * \li 2x2 complex "Jones" matrix
 * Note, there are two sets of Jones matrices, the second has terms multiplied 
 * in the reverse order of the first (multiplication doesn't generally commute).
 * \param in   Jones Calibrate Object.
 * \param cal  Calibrate Object.
 * \param sel  Data selector.
 * \param time desired time in days
 * \param SourID Desired Source ID
 * \param FreqID Desired Frequency ID
 * \param err  Error stack for messages and errors.
 */
static void ObitUVCalJonesUpdate (ObitUVCalJonesS *in, ObitUVCalCalibrateS *cal,
				  ObitUVSel *sel, ofloat time, olong SourID, olong FreqID,
				  ObitErr *err)
{
  olong i, iant, iif, ichan, index, jndex;
  gboolean good1, good2, bad;
  ofloat   wtt1, wtt2, delta, wt1, wt2;
  ofloat   wrk[8], fblank = ObitMagicF();
  ofloat   cpxblank[2] = {fblank,fblank};
  gboolean newcal;
  gchar *routine="ObitUVCalJonesUpdate";
 
 
  /* see if time for new table entry */
  if ((in->LastRowRead <= in->numRow)  &&  (time > in->FollowCalTime)) {
    ObitUVCalJonesNewTime (in, cal, sel, time, SourID, FreqID, err);
    if (err->error) Obit_traceback_msg (err, routine, "unspecified");
    newcal = TRUE;
  } else {
    newcal = FALSE;
  }

  /* see if calibration needs update; every 0.1 of solution interval. */
  delta = (time - in->CalTime);  
  if ((!newcal) &&  (delta <= 0.1*(in->FollowCalTime-in->PriorCalTime))) return;

  /* interpolate current calibration to time */
  in->CalTime    = time;
  in->curSourID  = SourID;
  /* initialize indices for Jones (index),  CalPrior and CalFollow (jndex) */
  index = 0;
  jndex = 0;

 /* loop thru antennas */
  for (iant=0; iant<in->numAnt; iant++) { /* antenna loop */

    /* set interpolation weights proportional to time difference. */
    wtt1 = 0.0;
    if (time < in->FollowAntTime[iant]) {
      if (in->FollowAntTime[iant] > in->PriorAntTime[iant]) 
	wtt1 = (in->FollowAntTime[iant] - time) / (in->FollowAntTime[iant] - in->PriorAntTime[iant]);
    } 
    wtt2 = 1.0 - wtt1;

    /* loop thru IF */
    for (iif=in->bIF; iif<=in->eIF; iif++) { /* IF loop */

      /* loop thru channels */
      for (ichan=in->bChan; ichan<=in->eChan; ichan++) { /* channel loop */

	/* First set of Jones matrices */
	/* initialize soln with blanks */
	ObitMatxSet(in->Jones[index], cpxblank, 0, 0);

	/* check for blanked soln. */
	good1 = (in->CalPrior[jndex] != fblank) && (in->CalPrior[jndex+1] != fblank) &&
	  (wtt1 > 0.0);
	good2 = (in->CalFollow[jndex] != fblank) &&  (in->CalFollow[jndex+1] != fblank)  && 
	  (wtt2 > 0.0);

	/* solution all flagged? */
	bad = !(good1 || good2);

	/* nothing more if both prior and following are bad */
	/* Interpolate Jones matrix, invert to in->Jones */
	if (!bad) {
	  
	  /* initial weights */
	  wt1 = wtt1; wt2 = wtt2;
	  
	  /* Only Following good */
	  if (!good1) { wt1 = 0.0; wt2 = 1.0; } 
	  
	  /* Only Prior good */
	  if (!good2) { wt1 = 1.0; wt2 = 0.0; } 

	  /* Average elements of Jones matrix, invert  */
	  for (i=0; i<8; i++) wrk[i] = wt1*in->CalPrior[jndex+i] + wt2*in->CalFollow[jndex+i];
	  ObitMatxSet2C(in->TempMatx, wrk[0],wrk[1],wrk[2],wrk[3],wrk[4],wrk[5],wrk[6],wrk[7]);
	  ObitMatxInv2x2(in->TempMatx, in->Jones[index]);
	} /* end of only valid solutions section */

	/* update indices */
        index += 1;  /* Array of ObitMatx inverse Jones matrix */
	jndex += 8;  /* Jones matrix as array of floats */

      } /* end channel loop */
    } /* end IF loop  */
  } /* end antenna loop */
} /* end ObitUVCalJonesUpdate */
    
/**
 * Read calibration for next time from JI/JT table.
 * If using the JI table, it is filled into both 
 CalPrior and CalFollow
 * Initialize on first call
 * \param in   Jones Calibrate Object.
 * \param cal  Calibrate Object.
 * \param sel  Data selector.
 * \param time desired time in days
 * \param SourID Desired Source ID
 * \param FreqID Desired Frequency ID
 * \param in   Error stack for messages and errors.
 */
 static void ObitUVCalJonesNewTime (ObitUVCalJonesS *in, ObitUVCalCalibrateS *cal, 
				    ObitUVSel *sel, ofloat time, olong SourID, olong FreqID,
				    ObitErr *err)
{
  ofloat fblank = ObitMagicF();
  olong nentry, i, j, ii, iant, iif, ichan, indx, jndx, lenEntry, lenEntryAnt;
  olong  ioff, joff, irow=0, limit, antno=1, nChIF;
  gboolean want=TRUE, done=FALSE;
  ObitTableJI *JITable = NULL;
  ObitTableJIRow *JITableRow = NULL;
  ObitTableJT *JTTable = NULL;
  ObitTableJTRow *JTTableRow = NULL;
  gchar *routine="ObitUVCalJonesNewTime";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;

  /* Sizes of things */
  nChIF = (in->eIF-in->bIF+1)*(in->eChan-in->bChan+1); /* Number of channels*IFs */
  lenEntry = 8;                    /* size of Jones matrix */
  lenEntryAnt = nChIF * lenEntry;  /* Size of CalPrior/CalFollow entry per antenna */
  
  /* initialize Prior and Following arrays if first call */
  if (in->LastRowRead <= 0) {
    /* Create */
    nentry = in->numAnt * (in->eIF-in->bIF+1) * (in->eChan-in->bChan+1) * lenEntry;
    in->CalPrior      = g_malloc0((nentry+3)*sizeof(ofloat));
    in->CalFollow     = g_malloc0((nentry+3)*sizeof(ofloat));
    in->PriorAntTime  = g_malloc0((in->numAnt+3)*sizeof(ofloat));
    in->FollowAntTime = g_malloc0((in->numAnt+3)*sizeof(ofloat));
    /* Blank fill */
    for (i=0; i<nentry; i++) in->CalPrior[i]  = fblank;
    for (i=0; i<nentry; i++) in->CalFollow[i] = fblank;
    for (i=0; i<in->numAnt; i++) in->PriorAntTime[i] = -1.0e10;
    for (i=0; i<in->numAnt; i++) in->FollowAntTime[i]= -1.0e10;
    /*for (i=0; i<in->numAnt; i++) cal->PriorIFR[i]     = fblank;
      for (i=0; i<in->numAnt; i++) cal->FollowIFR[i]    = fblank;*/
    in->PriorCalTime  = -1.0e10;
    in->FollowCalTime = -1.0e10;
    in->PriorSourID  = -1;
    in->FollowSourID = -1;
  } /* end of initialize on first call */

  /* Shuffle data from Following to Prior if time exceeded */
  for (iant=0; iant<in->numAnt; iant++) { /* antenna loop */
    /* if 2nd time exceeded - shift down */
    if ((time > in->FollowAntTime[iant])  &&  (in->FollowAntTime[iant] > -100.)) {
      in->PriorAntTime[iant] = in->FollowAntTime[iant];
      /*cal->PriorIFR[iant]    = cal->FollowIFR[iant]; NYI */
      for (iif=in->bIF; iif<=in->eIF; iif++) { /* IF loop */
	indx = lenEntryAnt * (iant) +  lenEntry * in->numChan * (iif-in->bIF);
	for (ichan=in->bChan; ichan<=in->eChan; ichan++) { /* channel loop */
	  jndx = indx + lenEntry * (ichan-in->bChan);
	  for (j=0; j<lenEntry; j++) in->CalPrior[jndx+j]  = in->CalFollow[jndx+j];
	} /* end channel loop  */
      } /* end IF loop  */
    } /* End if need to shuffle antenna */
  } /* end antenna loop */

  
  /* Jones table  - set local pointers */
  if (in->useJI) {
    JITable = (ObitTableJI*)in->JITable;
    JITableRow = (ObitTableJIRow*)in->JITableRow;
  } else if (in->useJT) {
    JTTable = (ObitTableJT*)in->JTTable;
    JTTableRow = (ObitTableJTRow*)in->JTTableRow;
  } else {
    Obit_log_error(err, OBIT_Error, "%s: Jones Matrix Table type not defined", routine);
    return;
  }
  
  /* Read through rows filling in data - read until selected time. */
  limit = MAX (1, in->LastRowRead);
  in->LastRowRead = 0;  /* The next time may start somewhere nonobvious */
  for (i=limit; i<=in->numRow; i++) { /* table row loop */
    irow = i;
    if (in->useJI) ObitTableJIReadRow (JITable, irow, JITableRow, err);
    if (in->useJT) ObitTableJTReadRow (JTTable, irow, JTTableRow, err);
    if (err->error) Obit_traceback_msg (err, routine, "Cal(Jones) table");
    if (in->useJI && (JITableRow->status < 0)) continue; /* JI entry flagged? */
    if (in->useJT && (JTTableRow->status < 0)) continue; /* JT entry flagged? */
    
    /* check subarray */
    if (in->useJI) {want = ((JITableRow->SubA == in->SubA)  ||  (JITableRow->SubA <= 0) || (in->SubA <= 0));}
    if (in->useJT) {want = ((JTTableRow->SubA == in->SubA)  ||  (JTTableRow->SubA <= 0) || (in->SubA <= 0));}
    
    /* check frqsel */
    if (in->useJI) {want = want &&
	((JITableRow->FreqID == FreqID) || (JITableRow->FreqID <= 0) || (FreqID <= 0));}
    if (in->useJT) {want = want &&
	((JTTableRow->FreqID == FreqID) || (JTTableRow->FreqID <= 0) || (FreqID <= 0));}
    
    /* check Source */
    if (in->useJI) {want = want && (ObitUVSelWantSour(sel, JITableRow->SourID) || (JITableRow->SourID==0));}
    if (in->useJT) {want = want && (ObitUVSelWantSour(sel, JTTableRow->SourID) || (JTTableRow->SourID==0));}
    
    /* Antenna number in range */
    if (in->useJI) {want = want && (JITableRow->antNo<=in->numAnt);}
    if (in->useJT) {want = want && (JTTableRow->antNo<=in->numAnt);}
    
    /* skip if not wanted */
    if (!want) continue;
    
    /* antenna number */
    if (in->useJI) antno = JITableRow->antNo;
    if (in->useJT) antno = JTTableRow->antNo;
    iant  = antno-1;
    
    /* time -> include this one? */
    if (time >= in->FollowAntTime[iant]) { 
      /* Copy Follow to Prior if empty */
      if (in->PriorAntTime[iant] > -100.) {
	/* new following entry - copy to prior */
	in->PriorAntTime[iant] = in->FollowAntTime[iant];
	in->PriorSourID        = in->FollowSourID;
	indx = lenEntryAnt * (iant);
	for (ii=0; ii<lenEntryAnt; ii++) in->CalPrior[indx+ii] = in->CalFollow[indx+ii];
      } /* end shuffle to Prior */
      
      /* fill in new following values */
      if (in->useJI) {
	in->FollowAntTime[iant] = JITableRow->Time;
	in->FollowCalTime       = JITableRow->Time;
	in->FollowSourID        = JITableRow->SourID;
      } else {
	in->FollowAntTime[iant] = JTTableRow->Time;
	in->FollowCalTime       = JTTableRow->Time;
	in->FollowSourID        = JTTableRow->SourID;
      } 
      /* copy entry from table row - select by channel,IF */
      for (iif=in->bIF; iif<=in->eIF; iif++) { /* IF loop */
	ioff = lenEntry * in->numChan * (iif-in->bIF);
	indx = lenEntryAnt*iant + lenEntry*(in->eChan-in->bChan+1)*(iif-in->bIF);
	for (ichan=in->bChan; ichan<=in->eChan; ichan++) { /* channel loop */
	  joff = ioff + lenEntry * (ichan-in->bChan);
	  jndx = indx + lenEntry * (ichan-in->bChan);
	  if (in->useJI) {
	    for (ii=0; ii<lenEntry; ii++) in->CalFollow[jndx+ii] = JITableRow->Jones[joff+ii];
	  }
	  if (in->useJT) {
	    for (ii=0; ii<lenEntry; ii++) in->CalFollow[jndx+ii] = JTTableRow->Jones[joff+ii];
	  }
	} /* end Channel loop */
      } /* end IF loop */      

      /* if Prior entry not valid copy following */
      if (in->PriorAntTime[iant] <= -100.) {
	in->PriorAntTime[iant] = in->FollowAntTime[iant];
	in->PriorSourID        = in->FollowSourID;
	indx = lenEntryAnt * (iant);
	for (ii=0; ii<lenEntryAnt; ii++) in->CalPrior[indx+ii]  = in->CalFollow[indx+ii];
      } /* end copy to Prior */
      
    } else {    
      /* This table entry not needed - are we there yet? */
      /* May need to restart earlier in table for some antennas.
	 Remember the first record not used. */
      
      if (in->LastRowRead <= 0) in->LastRowRead = i;
      
      /* if desired time for some antenna still after the Following time and there is no
	 Prior entry, keep going */
      done = TRUE;
      for (antno= 1; antno<=in->numAnt; antno++) { /* antennaloop */
	iant = antno-1;
	if ((time >= in->FollowAntTime[iant]) && (in->PriorAntTime[iant] >= -100.0)) done=FALSE;
      } /* end antenna loop */
      
      /* no more to fill in */
      if (done) break;
    } 
  } /* end loop over table entries */

  /* If all entries have been read and some antennas do not have a 
     following entry make dummy blanked entry at the last time */
  if (irow>=in->numRow ) {
    for (antno=1; antno<=in->numAnt; antno++) {
      iant = antno-1;
      if ((time >= in->FollowAntTime[iant]) && (in->PriorAntTime[iant] >= -100.0)) {
	/* Dummy Following entry */
	in->FollowAntTime[iant] += 10.0;
	/* Blank */
	indx = lenEntryAnt * (iant);
	for (ii=0; ii<lenEntryAnt; ii++) in->CalFollow[indx+ii] = fblank;
      } /* end if no entry */
    } /* end checking antennas  */
  } /* end if got to end */
  
  /* finished file using all entries? */
  if (in->LastRowRead <= 0) in->LastRowRead = in->numRow + 1;
  
  /* Set times */
  in->FollowCalTime = 1.0e10;
  in->PriorCalTime = -1.0e10;
  for (antno= 1; antno<=in->numAnt; antno++) { /* antenna loop*/
    iant = antno -1;
    if (in->PriorAntTime[iant] >= -100.0) {
      if (time >= in->PriorAntTime[iant]) 
	in->PriorCalTime = MAX (in->PriorCalTime, in->PriorAntTime[iant]);
      if (time <= in->FollowAntTime[iant]) 
	in->FollowCalTime = MIN (in->FollowCalTime,in->FollowAntTime[iant]);
    } 
  } /* end antenna loop */
  
  /* just to be sure something rational in times */
  if (in->PriorCalTime < -1000.0)  in->PriorCalTime  = time - 2.0/86400.0;
  if (in->FollowCalTime > 10000.0) in->FollowCalTime = time + 2.0/86400.0;
  return;
} /* end ObitUVCalJonesNewTime */


