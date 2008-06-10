/* $Id: ObitUVCalBaseline.c,v 1.3 2007/08/31 17:24:04 bcotton Exp $                            */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitUVCalBaseline.h"
#include "ObitUVCalBaselineDef.h"
#include "ObitUVDesc.h"
#include "ObitUVSel.h"
#include "ObitTableBL.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVCalBaseline.c
 * ObitUVCal utilities for applying baseline dependent calibration to
 * uv data.
 */

/*-------------------Private function prototypes-------------------*/
/** Private:  Create structure for baseline dependent Calibrate . */
static ObitUVCalBaselineS* newObitUVCalBaselineS (ObitUVCal *in);

/** Private: Update calibration arrays. */
static void ObitUVCalBaselineUpdate (ObitUVCalBaselineS *in, ofloat time,
					ObitErr *err);

/** Private:  Read calibration for a new time into the internal arrays. */
static void ObitUVCalBaselineNewTime (ObitUVCalBaselineS *in, ofloat time,
					ObitErr *err);

/*----------------------Public functions---------------------------*/
/**
 * Initialize structures for baseline dependent Calibrate .
 * \param in   Baseline Object.
 * \param sel  Data selector.
 * \param desc Data descriptor.
 * \param err  ObitError stack.
 */
void ObitUVCalBaselineInit (ObitUVCal *in, ObitUVSel *sel, ObitUVDesc *desc, 
		    ObitErr *err)
{
  ObitIOCode retCode;
  ObitUVCalBaselineS *me;
  olong size, i;
  gchar *routine="ObitUVCalBaselineInit";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVCalIsA(in));

  in->doBL = sel->doBLCal;
  if (!in->doBL) return;

  in->baselineCal = newObitUVCalBaselineS(in);

  /* pointer to calibration structure */
  me = in->baselineCal;

  /* Copy Selector information */
  me->bChan       = sel->startChann;
  me->eChan       = sel->startChann + sel->numberChann - 1;
  me->bIF         = sel->startIF;
  me->eIF         = sel->startIF + sel->numberIF - 1;
  me->SubA        = sel->SubA;
  me->FreqID      = sel->FreqID;

  /* Copy Cal information */
  me->BLTable     = ObitTableRef(in->BLTable);
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
    ObitTableBLOpen ((ObitTableBL*)(me->BLTable), OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);
  
  /* create row structure */
  me->BLTableRow = (Obit*)newObitTableBLRow((ObitTableBL*)(me->BLTable));
  
  /* table information */
  me->numPol = ((ObitTableBL*)me->BLTable)->numPol;
  me->numRow = ((ObitTableBL*)me->BLTable)->myDesc->nrow;

  /* Allocate calibration arrays */
  me->lenCalArrayEntry = 2; /* length of cal array entry */
  me->numBase = (me->numAnt)*(me->numAnt-1)/2; /* number of baselines */
  size = me->numBase * (me->eIF- me->bIF + 1) * me->numPol * me->lenCalArrayEntry;
  me->CalApply     = g_malloc(size*sizeof(float));
  me->CalPrior     = g_malloc(size*sizeof(float));
  me->CalFollow    = g_malloc(size*sizeof(float));

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


} /*  end ObitUVCalBaselineInit */

/**
 * Baseline data baseline dependent
 * A calibration array entry consists of:
 * \li real part of gain
 * \li imaginary of gain
 * \param in    Baseline Object.
 * \param time  Time of datum
 * \param ant1  first antenna number of baseline
 * \param ant2  second antanna of baseline.
 * \param RP    Random parameters array.
 * \param visIn 1 visibility as an array of floats
 * \param err   ObitError stack.
 */
void ObitUVCalBaseline (ObitUVCal *in, float time, olong ant1, olong ant2, 
			ofloat *RP, ofloat *visIn, ObitErr *err)
{
  olong   blndx, lentry, blindx, iif, ipol, ifreq, ioff, joff, index, nstoke, iSubA, itemp, numIF;
  gboolean   ccor;
  gboolean calBad;
  ofloat gwt, tvr, tvi, gr, gi, fblank = ObitMagicF();
  ObitUVCalBaselineS *me;
  ObitUVDesc *desc;
  ObitUVSel *sel;
  olong corID;
  gchar *routine="ObitUVCalBaseline";

  /* local pointers for structures */
  me   = in->baselineCal;
  desc = in->myDesc;
  sel  = in->mySel;

  /* Number of stokes correlations */
  nstoke = desc->inaxes[desc->jlocs];

  /* number of IFs */
  if (desc->jlocif>=0) numIF = desc->inaxes[desc->jlocif];
  else numIF = 1;

  /* Correlator ID if in data */
  if (desc->ilocid>=0) 
    corID = (olong)RP[desc->ilocid] + 0.5;
  else
    corID = 1;

  /* Subarray number in data */
  itemp = (olong)RP[desc->ilocb];
  iSubA = 1 + (olong)(100.0*(RP[desc->ilocb] -(ofloat)itemp) + 0.1);

  /* see if new time - update cal. */
  if ((time > me->CalTime) && (me->LastRowRead < me->numRow)) {
    ObitUVCalBaselineUpdate (me, time, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* check if cross correlation */
  ccor = ant1 != ant2;

  /* Set baseline index */
  blndx = ((ant1-1)*me->numAnt) - (((ant1+1)*ant1)/2) + ant2;
  lentry = 2 * (me->eIF - me->bIF + 1) * me-> numPol;
  blindx  = lentry * (blndx-1);

  /* loop over IF */
  for (iif= me->bIF; iif<=me->eIF; iif++) { /* loop 300 */
    ioff = (iif-1) * desc->incif;
    
    /* loop over polarization (only first two) */
    for (ipol= 1; ipol<=MAX(2,nstoke); ipol++) { /* loop 200 */

      /* check baseline flags */
      joff = ioff + (ipol-1) * desc->incs;

      /* Initialize calibration */
      gr  = 1.0;
      gi  = 0.0;
      gwt = 0.0;

      /* check IF flags */
      calBad = (me->CalApply[blindx] == fblank);

      /* Calibrate this one? */
      if (!calBad && ccor) {
	gr = me->CalApply[blindx];
	gi = me->CalApply[blindx+1];

 	/* "weight" calibration */
	if (me->doCalWt) {
	  gwt = (gr*gr + gi*gi);
	  if (gwt > 1.0e-10) gwt = 1.0 / gwt;
	} else {
	  gwt = 1.0;
	} 
      }

      /* loop over channels calibrating. */
      for (ifreq= me->bChan-1; ifreq<me->eChan; ifreq++) { /* loop 80 */
	index = joff + (ifreq) * desc->incf;
	tvr = gr * visIn[index]   + gi * visIn[index+1];
	tvi = gr * visIn[index+1] - gi * visIn[index];

	/* save calibrated data */
	visIn[index]   = tvr;
	visIn[index+1] = tvi;
	visIn[index+2] = visIn[index+2] * gwt;

      } /* end loop over channels  L80 */;

    } /* end loop loop over Stokes  L200 */;

  } /* end loop  L300 - loop over IF */;

} /* end ObitUVCalBaseline */


/**
 * Shutdown baseline dependent Calibrate.
 * Close any open file and destroy structures.
 * \param in   Baseline Object.
 * \param err  ObitError stack.
 */
void ObitUVCalBaselineShutdown (ObitUVCal *in, ObitErr *err)
{
  ObitUVCalBaselineS *me;
  ObitIOCode retCode;
  gchar *routine="ObitUVCalBaselineShutdown";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVCalIsA(in));

  /* Close calibration table  */
  me = in->baselineCal;
  retCode = ObitTableBLClose ((ObitTableBL*)me->BLTable, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);

  /* Release row structure  */
  me->BLTableRow = ObitTableBLRowUnref(me->BLTableRow);

  /* delete structure */
  in->baselineCal = ObitUVCalBaselineSUnref(in->baselineCal);
} /*  end ObitUVCalBaselineShutdown */

/**
 * Destroy structure for baseline dependent Calibrate .
 * \param in   Baseline Object.
 * \return NULL
 */
ObitUVCalBaselineS*
ObitUVCalBaselineSUnref (ObitUVCalBaselineS *in)
{
  if (in==NULL) return in;;
  in->BLTable    = ObitTableBLUnref((ObitTableBL*)in->BLTable);
  in->BLTableRow = ObitTableBLRowUnref((ObitTableBLRow*)in->BLTableRow);
  if (in->CalApply)     g_free(in->CalApply); in->CalApply   = NULL;
  if (in->CalPrior)     g_free(in->CalPrior); in->CalPrior   = NULL;
  if (in->CalFollow)    g_free(in->CalFollow); in->CalFollow  = NULL;

  /* basic structure */
   g_free (in);

  return NULL;
} /*  end ObitUVCalBaselineSUnref */

/*---------------Private functions---------------------------*/

/**
 * Create structure for baseline dependent Calibrate .
 * \param in   Baseline Object.
 * \return newly created object.
 */
static ObitUVCalBaselineS*
newObitUVCalBaselineS (ObitUVCal *in)
{
  ObitUVCalBaselineS* out;

  out = g_malloc0(sizeof(ObitUVCalBaselineS));

  /* Null pointers */
  out->BLTable    = NULL;
  out->BLTableRow = NULL;
  out->CalApply   = NULL;
  out->CalPrior   = NULL;
  out->CalFollow  = NULL;

  return out;
} /*  end newObitUVCalBaselineS */

/**
 * Update ObitUVCalBaseline calibration tables for time time.
 * The current table is interpolated between the previous and following
 * sets of solutions.
 * If a new set of entries is needed from the BL table they are read.
 * If the time for the FollowCalTime or PriorCalTime times are < -100 then 
 * it is considered that they do not contain valid data.
 * After the final BL table entries are past, or if there is only a single time,
 * the CalTime member is set to a large value so this routine will not be
 * called again.
 * Adopted from AIPS BLSET.FOR
 * A calibration array entry (prior or following) consists of:
 * \li real part of gain
 * \li imaginary of gain
 * \param in   Baseline Object.
 * \param time desired time in days
 * \param in   Error stace for messages and errors.
 */
static void ObitUVCalBaselineUpdate (ObitUVCalBaselineS *in, ofloat time,
					ObitErr *err)
{
  olong   loop, number;
  gboolean   good1, good2;
  ofloat      wtt1, wtt2, g1r, g1i, g2r, g2i, wt1, wt2, fblank = ObitMagicF();
  gchar *routine="ObitUVCalBaselineUpdate";
 
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;

  /* see if time for new table entry */
  if ((in->LastRowRead <= in->numRow)  &&  (time > in->FollowCalTime)) {
    ObitUVCalBaselineNewTime (in, time, err);
    if (err->error) Obit_traceback_msg (err, routine, "BL Table");
  }

  /* size of array */
  number = (in->numAnt * (in->numAnt-1)) / 2;
  number = number * (in->eIF - in->bIF + 1) * in->numPol * 2;
  
 /* set interpolation weights. */
  wtt1 = 0.0;
  if (in->FollowCalTime != in->PriorCalTime)
    wtt1 = (in->FollowCalTime - time) / (in->FollowCalTime - in->PriorCalTime);
  wtt2 = 1.0 - wtt1;

  /* Only Prior calibration? */
  if ((in->PriorCalTime>-100.0) && (in->FollowCalTime<-100.0)) {
    wtt1 = 1.0;
    wtt2 = 0.0;
  }
 
  /* Only Following calibration? */
  if ((in->PriorCalTime<-100.0) && (in->FollowCalTime>-100.0)) {
    wtt1 = 0.0;
    wtt2 = 1.0;
  }
 
  for (loop= 0; loop<number; loop += in->lenCalArrayEntry) { /* loop 500 */
    /* initialize soln with blanks */
    in->CalApply[loop]   = fblank;
    in->CalApply[loop+1] = fblank;
    
    /* check for blanked soln. */
    good1 = (in->CalPrior[loop] != fblank)  &&  (in->CalPrior[loop+1] != fblank)  && 
      (wtt1 > 0.0);
    good2 = (in->CalFollow[loop] != fblank)  &&  (in->CalFollow[loop+1] != fblank)  && 
      (wtt2 > 0.0);

    /* anything usable? */
    if (good1 || good2) {
      wt1 = wtt1;
      wt2 = wtt2;
      if (!good1) wt1 = 0.0;
      if (!good1) wt2 = 1.0;
      if (!good2) wt1 = 1.0;
      if (!good2) wt2 = 0.0;
 
     /* set values */
      g1r = in->CalPrior[loop];
      g1i = in->CalPrior[loop+1];
      g2r = in->CalFollow[loop];
      g2i = in->CalFollow[loop+1];
      
      /* interpolate */
      in->CalApply[loop]   = wt1 * g1r + wt2 * g2r;
      in->CalApply[loop+1] = wt1 * g1i + wt2 * g2i;
    }
  } /* end loop  L500: */;
  
  /* If time is still past Following and all the BL table has been read, 
     then there is nothing more to gain;  set following time to large number
     so this routine not called again. */
  if ((time>in->FollowCalTime) && ((in->LastRowRead >= in->numRow)))
      in->CalTime = 1.0e20;

} /* end ObitUVCalBaselineUpdate */
    
/**
 * Read calibration for next time from BL table.
 * If the time for the FollowCalTime or PriorCalTime times are < -100 then 
 * it is considered that they do not contain valid data.
 * Adopted from AIPS BLGET.FOR
 * \param in   Calibrate Object.
 * \param time desired time in days
 * \param in   Error stack for messages and errors.
 */
static void ObitUVCalBaselineNewTime (ObitUVCalBaselineS *in, ofloat time,
					ObitErr *err)
{
  ObitIOCode retCode;
  ObitTableBL *BLTable = NULL;
  ObitTableBLRow *BLTableRow = NULL;
  ofloat *ftemp, ft, curtim=-1.0, fblank = ObitMagicF();
  olong nblank, i, irow, limit;
  gboolean want, got1;
  olong   ifno, lentry, ant1, ant2, blndx, ifoff;
  olong   index,  numbl;
  gchar *routine="ObitUVCalBaselineNewTime";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;

  /* Have we already read all of BL Table */
  if (in->LastRowRead >= in->numRow) return;

  /* BL table - set local pointers */
  BLTable = (ObitTableBL*)in->BLTable;
  BLTableRow = (ObitTableBLRow*)in->BLTableRow;

  /* setup on first call */
  if (in->LastRowRead <= 0) {

    /* blank fill arrays */
    numbl = (in->numAnt * (in->numAnt-1)) / 2;
    nblank = numbl * (in->eIF - in->bIF + 1) * in->numPol * in->lenCalArrayEntry;
    for (i= 0; i<nblank; i++) {
      in->CalPrior[i]  = fblank;
      in->CalFollow[i] = fblank;
    } 
    
    /* read until selected time. */
    lentry = in->lenCalArrayEntry * in->numPol * (in->eIF - in->bIF + 1);
    got1 = FALSE;
    limit = MAX (1, in->LastRowRead+1);

    /* loop over table */
    for (i= limit; i<=in->numRow; i++) { /* loop 80 */
      irow = i;
      retCode = ObitTableBLReadRow (BLTable, irow, BLTableRow, err);
      in->LastRowRead = i;
      if (err->error) Obit_traceback_msg (err, routine, "BL Table");
      if (BLTableRow->status < 0) continue; /* entry flagged? */
 
      /* check subarray */
      want = (BLTableRow->SubA == in->SubA) || (BLTableRow->SubA <= 0) || (in->SubA <= 0);
      
      /* check frqsel */
      want = want &&
	((BLTableRow->FreqID != in->FreqID) || (BLTableRow->FreqID <= 0) || (in->FreqID >= 0));
      
      if (!want) break; /* want this one? */
      
      /* see if time past data time. */
      if (BLTableRow->Time > time) break;
      
      /* fill in values */
      got1 = TRUE;
      in->PriorCalTime = BLTableRow->Time;
      
      /* which antennas */
      ant1 = BLTableRow->ant1;
      ant2 = BLTableRow->ant2;

      blndx = ((ant1-1)*in->numAnt) - ((ant1+1)*ant1)/2 + ant2;
      index = lentry * (blndx - 1);

     /* loop over IF */
      for (ifno= in->bIF; ifno<=in->eIF; ifno++) { /* loop 70 */
	ifoff = ifno - 1;
	in->CalPrior[index]   = BLTableRow->RealM1[ifoff];
	in->CalPrior[index+1] = BLTableRow->ImagM1[ifoff];
	index += in->lenCalArrayEntry;

	if (in->numPol > 1) {
	  /* second polarization */
	  in->CalPrior[index]   = BLTableRow->RealM2[ifoff];
	  in->CalPrior[index+1] = BLTableRow->ImagM2[ifoff];
	  index += in->lenCalArrayEntry;
	} 
      } /* end IF loop  L70:  */;
    } /* end loop  over table rows L80:  */
    /* end of initial reading */
  } else {

    /* not first call, swap Prior and Following */
    ftemp = in->CalPrior;
    in->CalPrior = in->CalFollow;
    in->CalFollow = ftemp;
    ft = in->PriorCalTime;
    in->PriorCalTime = in->FollowCalTime;
    in->FollowCalTime = ft;
  }

  /* Look for entries following time */

  /* blank fill CalFollow */
  numbl = (in->numAnt * (in->numAnt-1)) / 2;
  nblank = numbl * (in->eIF - in->bIF + 1) * in->numPol * in->lenCalArrayEntry;
  for (i= 0; i<nblank; i++) in->CalFollow[i] = fblank;
  
  /* check if bl table exhausted, if so just mark tables and quit */
  if (in->LastRowRead >= in->numRow) {
    in->FollowCalTime = -1.0e20;
    return;
  }

  /* read until selected time. */
  lentry = in->lenCalArrayEntry * in->numPol * (in->eIF - in->bIF + 1);
  limit = MAX (1, in->LastRowRead+1);

    /* loop over table */
    for (i= limit; i<=in->numRow; i++) { /* loop 180 */
      irow = i;
      retCode = ObitTableBLReadRow (BLTable, irow, BLTableRow, err);
      in->LastRowRead = i;
      if (err->error) Obit_traceback_msg (err, routine, "BL Table");
      if (BLTableRow->status < 0) continue; /* entry flagged? */
 
      /* check subarray */
      want = (BLTableRow->SubA == in->SubA)  ||  (BLTableRow->SubA <= 0) || (in->SubA <= 0);
      
      /* check frqsel */
      want = want &&
	((BLTableRow->FreqID != in->FreqID)  ||  (BLTableRow->FreqID <= 0)  || (in->FreqID >= 0));
      
      if (!want) break; /* want this one? */

      /* Want just the first time found */
      if (i == limit) curtim = BLTableRow->Time;
      
      /* see if time past desired time. */
      if (BLTableRow->Time > curtim) {
	in->LastRowRead--;   /* Start with this next time */
	break;
      }
      
      /* fill in values */
      in->FollowCalTime = BLTableRow->Time;
      
      /* which antennas */
      ant1 = BLTableRow->ant1;
      ant2 = BLTableRow->ant2;

      blndx = ((ant1-1)*in->numAnt) - ((ant1+1)*ant1)/2 + ant2;
      index = lentry * (blndx - 1);

     /* loop over IF */
      for (ifno= in->bIF; ifno<=in->eIF; ifno++) { /* loop 170 */
	ifoff = ifno - 1;
	in->CalFollow[index]   = BLTableRow->RealM1[ifoff];
	in->CalFollow[index+1] = BLTableRow->ImagM1[ifoff];
	index += in->lenCalArrayEntry;

	if (in->numPol > 1) {
	  /* second polarization */
	  in->CalFollow[index]   = BLTableRow->RealM2[ifoff];
	  in->CalFollow[index+1] = BLTableRow->ImagM2[ifoff];
	  index += in->lenCalArrayEntry;
	} 
      } /* end IF loop  L170:  */;
    } /* end loop  over table rows L180:  */

} /* end ObitUVCalBaselineNewTime */

