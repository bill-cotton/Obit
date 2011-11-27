/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2011                                          */
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

#include <glib.h>
#include "ObitUVCalFlag.h"
#include "ObitUVCalFlagDef.h"
#include "ObitUVDesc.h"
#include "ObitUVSel.h"
#include "ObitTableFG.h"
#include "ObitTableUtil.h"
#include "ObitThread.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVCalFlag.c
 * ObitUVCal utilities for applying flagging to uv data 
 */

/*-------------------Private function prototypes-------------------*/
/** Private:  Create structure for flagging. */
static ObitUVCalFlagS* newObitUVCalFlagS (ObitUVCal *in);

/** Private: Update flagging arrays. */
static void ObitUVCalFlagUpdate (ObitUVCalFlagS *in, ObitUVSel *sel, ofloat time,
				 ObitErr *err);

/** Private: Make Threaded args */
static olong MakeCalFlagFuncArgs (ObitThread *thread, ObitUVCalFlagS *in,
				  ObitUVDesc *desc,
				  CalFlagFuncArg ***ThreadArgs);

/** Private: Delete Threaded args */
static void KillCalFlagFuncArgs (olong nargs, CalFlagFuncArg **ThreadArgs);

/** Private: Threaded Flag routine */
static gpointer ThreadCalFlag (gpointer arg);

/** Private: Flag everything in a vis */
static void CalFlagAll (ObitUVDesc *desc, ofloat *visIn);

/** Private: Time to String */
static void T2String (ofloat time, gchar *msgBuf);
/*----------------------Public functions---------------------------*/
/**
 * Initialize structures for Flagging.
 * \param in   Flag Object.
 * \param sel  Data selector.
 * \param desc Data descriptor.
 * \param err  ObitError stack.
 */
void ObitUVCalFlagInit (ObitUVCal *in, ObitUVSel *sel, ObitUVDesc *desc, 
			ObitErr *err)
{
  ObitIOCode retCode;
  ObitUVCalFlagS *me;
  gchar *colName="TIME RANGE";
  gchar *routine="ObitUVCalFlagInit";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVCalIsA(in));

  in->doFlag = sel->doFlag;
  if (!in->doFlag) return;

  in->flag = newObitUVCalFlagS(in);

  /* pointer to flagging structure */
  me = in->flag;

  /* Copy Selector information */
  me->SubA        = sel->SubA;
  me->FreqID      = sel->FreqID;

  /* Copy Cal information */
  me->FGTable     = ObitTableRef(in->FGTable);
  me->LastRowRead = 0;
  me->numSubA     = in->numSubA;
  me->numIF       = desc->inaxes[desc->jlocif];
  me->numChan     = desc->inaxes[desc->jlocf];
  me->numStok     = desc->inaxes[desc->jlocs];
 
 /* First stokes parameter in data */
  if (desc->crval[desc->jlocs] > 0) 
    me->stoke0 = desc->crval[desc->jlocs] + 
      (1.0-desc->crpix[desc->jlocs]) * desc->cdelt[desc->jlocs]  + 0.5;
  else
    me->stoke0 = desc->crval[desc->jlocs] + 
      (1.0-desc->crpix[desc->jlocs]) * desc->cdelt[desc->jlocs]  - 0.5;

  /* Copy descriptor information */
  me->numAnt    = desc->maxAnt;
  me->numSubA   = desc->numSubA;
  me->numIF     = desc->inaxes[desc->jlocif];
  me->numChan   = desc->inaxes[desc->jlocf];

  /* Open Flagging table  */
  /* Sort to time order if needed */
  retCode = ObitTableUtilSort((ObitTable*)(me->FGTable), colName, FALSE, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);
  retCode = 
    ObitTableFGOpen ((ObitTableFG*)(me->FGTable), OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);
  
  /* create row structure */
  me->FGTableRow = (Obit*)newObitTableFGRow((ObitTableFG*)(me->FGTable));
  
  /* table information */
  me->numRow = ((ObitTableFG*)me->FGTable)->myDesc->nrow;

  /* Allocate flagging arrays */
  me->maxFlag     = 250000;
  me->numFlag     = 0;
  me->flagSour    = g_malloc0(me->maxFlag*sizeof(olong));
  me->flagAnt     = g_malloc0(me->maxFlag*sizeof(olong));
  me->flagBase    = g_malloc0(me->maxFlag*sizeof(olong));
  me->flagSubA    = g_malloc0(me->maxFlag*sizeof(olong));
  me->flagFQID    = g_malloc0(me->maxFlag*sizeof(olong));
  me->flagBIF     = g_malloc0(me->maxFlag*sizeof(olong));
  me->flagEIF     = g_malloc0(me->maxFlag*sizeof(olong));
  me->flagBChan   = g_malloc0(me->maxFlag*sizeof(olong));
  me->flagEChan   = g_malloc0(me->maxFlag*sizeof(olong));
  me->flagPol     = g_malloc0(4*me->maxFlag*sizeof(gboolean));
  me->flagEndTime = g_malloc0(me->maxFlag*sizeof(ofloat));

  /* Initial time to trigger update of calibration arrays */
  me->flagTime   = -1.0e20;
  me->maxSimFlag = 0;  /* Max. simultaneous flags */
  me->flagAll = FALSE;
  /* Setup threading structures */
  me->nThArg = MakeCalFlagFuncArgs (in->thread, me, desc, 
				    (CalFlagFuncArg***)&me->thArgArr);


} /*  end ObitUVCalFlagInit */

/**
 * Flag a visibility, possibly with threading
 * Adapted from AIPS DATFLG.FOR
 * \param in    Flag Object.
 * \param time  Time of datum
 * \param ant1  first antenna number of baseline
 * \param ant2  second antanna of baseline.
 * \param RP    Random parameters array.
 * \param visIn visibility as an array of floats
 * \param err   ObitError stack.
 */
void ObitUVCalFlag (ObitUVCal *in, float time, olong ant1, olong ant2, 
		    ofloat *RP, ofloat *visIn, ObitErr *err)
{
  olong   kbase, FQID, SourID, iSubA;
  olong   i, nElem, nElemPerThread, nTh, loElem, hiElem;
  gboolean OK;
  ObitUVCalFlagS *me;
  ObitUVDesc *desc;
  ObitUVSel *sel;
  CalFlagFuncArg **threadArgs;
  gchar *routine="ObitUVCalFlag";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVCalIsA(in));

  /* local pointers for structures */
  me   = in->flag;
  desc = in->myDesc;
  sel  = in->mySel;

  /* see if new time - update cal. */
  if (time > me->flagTime) {
    ObitUVCalFlagUpdate (me, sel, time, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* If no currently active flags, nothing to do */
  if (me->numFlag <= 0) return;
  
  /* Baseline and subarray number in data */
  kbase = (olong)RP[desc->ilocb];
  iSubA = 1 + (olong)(100.0*(RP[desc->ilocb] -(ofloat)kbase) + 0.1);

  /* Data FQ id */
  if (desc->ilocfq >= 0) FQID = RP[desc->ilocfq] + 0.1;
  else  FQID = 0;

  /* Source ID */
  if (desc->ilocsu >= 0) SourID = RP[desc->ilocsu] + 0.1;
  else SourID = 0;

  /* Are we flagging everything? */
  if (me->flagAll) {
    CalFlagAll (desc, visIn);
    return;
  }

  /* Divide up work, can have me->nThArg threads */
  nElem = me->numFlag;
  nElemPerThread = nElem/me->nThArg;
  nTh = me->nThArg;
  /* Lower limit on threading */
  if (nElem<1000) {nElemPerThread = nElem; nTh = 1;}
  loElem = 1;
  hiElem = nElemPerThread;
  hiElem = MIN (hiElem, nElem);

  /* Set up thread arguments */
  threadArgs = (CalFlagFuncArg**)me->thArgArr;
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hiElem = nElem;  /* Make sure do all */
    threadArgs[i]->first   = loElem-1;
    threadArgs[i]->last    = hiElem-1;
    threadArgs[i]->iSubA   = iSubA;
    threadArgs[i]->kbase   = kbase;
    threadArgs[i]->SourID  = SourID;
    threadArgs[i]->FQID    = FQID;
    threadArgs[i]->ant1    = ant1;
    threadArgs[i]->ant2    = ant2;
    threadArgs[i]->visIn   = visIn;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which Elems */
    loElem += nElemPerThread;
    hiElem += nElemPerThread;
    hiElem = MIN (hiElem, nElem);
  }

  /* Do operation */
  OK = ObitThreadIterator (in->thread, nTh, 
			   (ObitThreadFunc)ThreadCalFlag,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) {
    Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
  }
} /* end ObitUVCalFlag */


/**
 * Shutdown baseline dependent Calibrate.
 * Close any open file and destroy structures.
 * \param in   Flag Object.
 * \param err  ObitError stack.
 */
void ObitUVCalFlagShutdown (ObitUVCal *in, ObitErr *err)
{
  ObitUVCalFlagS *me;
  ObitIOCode retCode;
  gchar *routine="ObitUVCalFlagShutdown";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVCalIsA(in));

  /* Close calibration table  */
  me = in->flag;
  retCode = ObitTableFGClose ((ObitTableFG*)me->FGTable, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);

  /* Release row structure  */
  me->FGTableRow = ObitTableFGRowUnref(me->FGTableRow);

  /* If err->prtLv>2 tell max number of simultaneous flags */
  if ((err->prtLv>2) && (me->maxSimFlag>0)) {
    Obit_log_error(err, OBIT_InfoErr, "Maximum simultaneous flags %d", 
		   me->maxSimFlag);
  }

  /* delete structure */
  in->flag = ObitUVCalFlagSUnref(in->flag);
} /*  end ObitUVCalFlagShutdown */

/**
 * Destroy structure for baseline dependent Calibrate .
 * \param in   Flag Object.
 * \return NULL
 */
ObitUVCalFlagS*
ObitUVCalFlagSUnref (ObitUVCalFlagS *in)
{
  if (in==NULL) return in;;
  in->FGTable    = ObitTableFGUnref((ObitTableFG*)in->FGTable);
  in->FGTableRow = ObitTableFGRowUnref((ObitTableFGRow*)in->FGTableRow);
  if (in->flagSour)    g_free(in->flagSour);    in->flagSour  = NULL;
  if (in->flagAnt)     g_free(in->flagAnt);     in->flagAnt   = NULL;
  if (in->flagBase)    g_free(in->flagBase);    in->flagBase  = NULL;
  if (in->flagSubA)    g_free(in->flagSubA);    in->flagSubA  = NULL;
  if (in->flagFQID)    g_free(in->flagFQID);    in->flagFQID  = NULL;
  if (in->flagBIF)     g_free(in->flagBIF);     in->flagBIF   = NULL;
  if (in->flagEIF)     g_free(in->flagEIF);     in->flagEIF   = NULL;
  if (in->flagBChan)   g_free(in->flagBChan);   in->flagBChan = NULL;
  if (in->flagEChan)   g_free(in->flagEChan);   in->flagEChan = NULL;
  if (in->flagPol)     g_free(in->flagPol);     in->flagPol   = NULL;
  if (in->flagEndTime) g_free(in->flagEndTime); in->flagEndTime = NULL;
  if (in->thArgArr) {
    KillCalFlagFuncArgs (in->nThArg, in->thArgArr);
    in->thArgArr = NULL;
    in->nThArg = 0;
  }

  /* basic structure */
   g_free (in);

  return NULL;
} /*  end ObitUVCalFlagSUnref */

/*---------------Private functions---------------------------*/

/**
 * Create structure for baseline dependent Calibrate .
 * \param in   Flag Object.
 * \return newly created object.
 */
static ObitUVCalFlagS*
newObitUVCalFlagS (ObitUVCal *in)
{
  ObitUVCalFlagS* out;

  out = g_malloc0(sizeof(ObitUVCalFlagS));

  /* Null pointers */
  out->FGTable     = NULL;
  out->FGTableRow  = NULL;
  out->flagSour    = NULL;
  out->flagAnt     = NULL;
  out->flagBase    = NULL;
  out->flagSubA    = NULL;
  out->flagFQID    = NULL;
  out->flagBIF     = NULL;
  out->flagEIF     = NULL;
  out->flagBChan   = NULL;
  out->flagEChan   = NULL;
  out->flagPol     = NULL;
  out->flagEndTime = NULL;
  out->thArgArr    = NULL;
  
  return out;
} /*  end newObitUVCalFlagS */

/**
 * Update ObitUVCalFlag calibration tables for time time.
 * Adapted from AIPS NXTFLG.FOR
 * \param in   Flag Object.
 * \param in   UV selector object
 * \param time desired time in days
 * \param in   Error stack for messages and errors.
 */
static void ObitUVCalFlagUpdate (ObitUVCalFlagS *in, ObitUVSel *sel, ofloat time,
				 ObitErr *err)
{
  olong   ndrop, mdrop, a1, a2, it;
  olong irow,  i, limit4, nwords;
  gboolean done, dropall, warn;
  ObitIOCode retCode;
  ObitTableFG *FGTable = NULL;
  ObitTableFGRow *FGTableRow = NULL;
  gchar tString[25];
  size_t nbytes;
  gchar *routine="ObitUVCalFlagUpdate";
 
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert(in!=NULL);
  g_assert(ObitUVSelIsA(sel));

  /* New time for flags */
  in->flagTime = time;
 
  /* check if any flags expired. */
  done = (in->numFlag <= 0);
  while (!done) {
  
    /* find highest number expired flag */
    ndrop = -1;
    for (i=in->numFlag-1; i>=0; i--) { /* loop 100 */
      if (in->flagEndTime[i] < time) {ndrop = i; break;}
    } /* end loop  L100: */

    /* see if any to be dropped. */
    done = (ndrop < 0);
    if (done) break;

    /* See if any lower consecutive entries to be dropped */
    mdrop = 1; 
    i = ndrop;
    while (((i-mdrop)>=0) && (in->flagEndTime[i-mdrop] < time)) {
      mdrop++;
      ndrop--;
    }


    /* Are all dropped to end of list? */
    dropall = ((ndrop+mdrop) >= in->numFlag);
    if (dropall) {in->numFlag -= mdrop; continue;}

    /* compress, dropping flag. */
    if ((ndrop+mdrop) < in->numFlag) {  /* This is not the last one the list */
      i      = ndrop;
      nwords = in->numFlag - ndrop - mdrop + 1;
      nbytes = nwords * sizeof(ofloat);
      memmove(&in->flagEndTime[i], &in->flagEndTime[i+mdrop], nbytes);
      nbytes = nwords * sizeof(olong);
      memmove(&in->flagSour[i] , &in->flagSour[i+mdrop],  nbytes);
      memmove(&in->flagAnt[i]  , &in->flagAnt[i+mdrop],   nbytes);
      memmove(&in->flagFQID[i] , &in->flagFQID[i+mdrop],  nbytes);
      memmove(&in->flagBase[i] , &in->flagBase[i+mdrop],  nbytes);
      memmove(&in->flagSubA[i] , &in->flagSubA[i+mdrop],  nbytes);
      memmove(&in->flagBIF[i]  , &in->flagBIF[i+mdrop],   nbytes);
      memmove(&in->flagEIF[i]  , &in->flagEIF[i+mdrop],   nbytes);
      memmove(&in->flagBChan[i], &in->flagBChan[i+mdrop], nbytes);
      memmove(&in->flagEChan[i], &in->flagEChan[i+mdrop], nbytes);
      nbytes = 4 * nwords * sizeof(gboolean);
      memmove(&in->flagPol[i*4], &in->flagPol[4*(i+mdrop)], nbytes);

      /*    for (i= limit; i<in->numFlag-1; i++) { 
	    in->flagEndTime[i] = in->flagEndTime[i+1];
	    in->flagSour[i]    = in->flagSour[i+1];
	    in->flagAnt[i]     = in->flagAnt[i+1];
	    in->flagFQID[i]    = in->flagFQID[i+1];
	    in->flagBase[i]    = in->flagBase[i+1];
	    in->flagSubA[i]    = in->flagSubA[i+1];
	    in->flagBIF[i]     = in->flagBIF[i+1];
	    in->flagEIF[i]     = in->flagEIF[i+1];
	    in->flagBChan[i]   = in->flagBChan[i+1];
	    in->flagEChan[i]   = in->flagEChan[i+1];
	    in->flagPol[i*4]   = in->flagPol[4*(i+1)];
	    in->flagPol[i*4+1] = in->flagPol[4*(i+1)+1];
	    in->flagPol[i*4+2] = in->flagPol[4*(i+1)+2];
	    in->flagPol[i*4+3] = in->flagPol[4*(i+1)+3];
	    }*/
    }  
    in->numFlag -= mdrop;
  } /* end loop deleting expired entries */

    
  /* check if list exhausted */
  if (in->LastRowRead >= in->numRow) goto alldone;

  /* FG table - set local pointers */
  FGTable = (ObitTableFG*)in->FGTable;
  FGTableRow = (ObitTableFGRow*)in->FGTableRow;

  warn = TRUE;
  /* Find next valid flag. Loop through records */
  limit4 = MAX (1, in->LastRowRead+1);
  for (i= limit4; i<=in->numRow; i++) { /* loop 360 */
    irow = i;
    retCode = ObitTableFGReadRow (FGTable, irow, FGTableRow, err);
    if (err->error) Obit_traceback_msg (err, routine, "FG Table");
    if (FGTableRow->status < 0) continue; /* entry flagged? */
    in->LastRowRead = i-1;  /* may try this one again */
  
    if (time < FGTableRow->TimeRange[0]) goto alldone;
    if (time > FGTableRow->TimeRange[1]) continue;
  
    /* check FQ id. */
    if ((FGTableRow->freqID != in->FreqID)  &&  (in->FreqID > 0)  && 
	(FGTableRow->freqID > 0)) continue;

    /* Check that starting IF is in range */
    if ((FGTableRow->ifs[0] > 0) &&  	(FGTableRow->ifs[0] > in->numIF)) continue;

    /* Check that starting channel is in range */
    if ((FGTableRow->chans[0] > 0) && (FGTableRow->chans[0] > in->numChan)) continue;
    
    /* Is this source wanted */
    if ((FGTableRow->SourID>0) && 
	(!ObitUVSelWantSour(sel, (olong)FGTableRow->SourID))) continue;
    
    /* Must want this one - save it */
    /* check if too big */
    if (in->numFlag >= in->maxFlag) {
      if (warn) {
	warn = FALSE; /* Only once */
	Obit_log_error(err, OBIT_InfoWarn, 
		       "ERROR: Exceeded limit of %d simultaneous flags for %s", 
		       in->maxFlag, sel->name);
	T2String (time, tString);
	Obit_log_error(err, OBIT_InfoWarn, 
		       "       Time: %s", tString);
	Obit_log_error(err, OBIT_InfoWarn, 
		       "       Flagging all data until number drops");
	ObitErrLog(err); 
	in->flagAll = (in->numFlag >= in->maxFlag);
      }
      continue;
    }

    /* New flag */
    in->numFlag++;
 
    /* fill in tables */
    in->LastRowRead = i;  /* Last row actually used */
    it = in->numFlag - 1;
    in->flagEndTime[it] = FGTableRow->TimeRange[1];
    in->flagSour[it] = FGTableRow->SourID;
    in->flagFQID[it] = FGTableRow->freqID;
    a1 = MIN (FGTableRow->ants[0], FGTableRow->ants[1]);
    a2 = MAX (FGTableRow->ants[0], FGTableRow->ants[1]);
    if (a1 <= 0) {
      in->flagAnt[it] = a2;
      in->flagBase[it] = 0;
    } else {
      in->flagAnt[it]  = a1;
      in->flagBase[it] = a1*256 + a2;
    } 
    in->flagSubA[it] = FGTableRow->SubA;
    in->flagBIF[it] = FGTableRow->ifs[0];
    in->flagEIF[it] = FGTableRow->ifs[1];
    if (in->flagBIF[it] <= 0) in->flagBIF[it] = 1;
    if (in->flagEIF[it] <= 0) in->flagEIF[it] = in->numIF;

    in->flagBChan[it] = FGTableRow->chans[0];
    in->flagEChan[it] = MIN (in->numChan, FGTableRow->chans[1]);
    if (in->flagBChan[it] <= 0) in->flagBChan[it] = 1;
    if (in->flagEChan[it] <= 0) in->flagEChan[it] = in->numChan;

    /* Ensure that IF and channel selection are in range */
    in->flagBIF[it]   = MAX (in->flagBIF[it], 1);
    in->flagBChan[it] = MAX (in->flagBChan[it], 1);
    in->flagEIF[it]   = MIN (in->flagEIF[it], in->numIF);
    in->flagEChan[it] = MIN (in->flagEChan[it], in->numChan);

    /* Stokes flags are packed into a bit array */
    in->flagPol[it*4]   = FGTableRow->pFlags[0] & 0x1;
    in->flagPol[it*4+1] = FGTableRow->pFlags[0] & 0x2;
    in->flagPol[it*4+2] = FGTableRow->pFlags[0] & 0x4;
    in->flagPol[it*4+3] = FGTableRow->pFlags[0] & 0x8;
    
  }  /* end loop over file  L360: */

 alldone:
  /* Save maximum number of simultaneous flags */
  in->maxSimFlag = MAX (in->maxSimFlag, in->numFlag);
  /* Resuming normal flagging */
  if (in->flagAll && ((in->numFlag<in->maxFlag))) {
    T2String (time, tString);
    Obit_log_error(err, OBIT_InfoWarn, 
		   "       Resume normal flagging at Time: %s", tString);
    ObitErrLog(err); 
  }
  in->flagAll = (in->numFlag >= in->maxFlag);

} /* end ObitUVCalFlagUpdate */

/**
 * Make arguments for a Threaded ThreadCalFlag
 * \param thread     ObitThread object to be used
 * \param in         ObitUVCalFlagS to be operated on
 * \param desc       ObitUVDesc for data
 * \param ThreadArgs[out] Created array of FAFuncArg, 
 *                   delete with KillFAFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
static olong MakeCalFlagFuncArgs (ObitThread *thread, ObitUVCalFlagS *in,
				  ObitUVDesc *desc,
				  CalFlagFuncArg ***ThreadArgs)
{
  olong i, nThreads;

  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(thread));

  /* Initialize threadArg array */
  *ThreadArgs = g_malloc0(nThreads*sizeof(CalFlagFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*ThreadArgs)[i] = g_malloc0(sizeof(CalFlagFuncArg)); 
  for (i=0; i<nThreads; i++) {
    (*ThreadArgs)[i]->thread = ObitThreadRef(thread);
    (*ThreadArgs)[i]->in     = in;
    (*ThreadArgs)[i]->desc   = ObitUVDescRef(desc);
    (*ThreadArgs)[i]->first = 0;
    (*ThreadArgs)[i]->last  = 0;
    (*ThreadArgs)[i]->ithread  = i;
  }

  return nThreads;
} /*  end MakeCalFlagFuncArg */

/**
 * Delete arguments for ThreadCalFlag
 * \param nargs      number of elements in ThreadArgs.
 * \param ThreadArgs Array of CalFlagFuncArg
 */
static void KillCalFlagFuncArgs (olong nargs, 
				 CalFlagFuncArg **ThreadArgs)
{
  olong i;

  if (ThreadArgs==NULL) return;
  ObitThreadPoolFree (ThreadArgs[0]->thread);  /* Free thread pool */
  for (i=0; i<nargs; i++) {
    if (ThreadArgs[i]) {
      if (ThreadArgs[i]->thread) ObitThreadUnref(ThreadArgs[i]->thread);
      if (ThreadArgs[i]->desc)   ObitUVDescUnref(ThreadArgs[i]->desc);
      g_free(ThreadArgs[i]);
    }
  }
  g_free(ThreadArgs);
} /*  end KillCalFlagFuncArgs */

/**
 * Thread apply flags to data
 * Callable as thread
 * The potential dependencies between threads should not be a problem.
 * \param arg Pointer to CalFlagFuncArg argument with elements:
 * \li   thread ObitThread to use 
 * \li   in     CalFlag to work on 
 * \li   first  First flag entry (0-rel) number 
 * \li   last   Highest flag entry  (0-rel) number 
 * \li   time   Time of datum (day) 
 * \li   iSubA  Subarray number  of datum
 * \li   kbase  Baseline index  of datum 
 * \li   SourID Source ID of datum
 * \li   FQID   FQ ID of datum
 * \li   ant1   First antenna number
 * \li   ant2   Second antenna number
 * \li   visIn  Visibility as an array of floats 
 * \li   ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadCalFlag (gpointer arg)
{
  /* Get arguments from structure */
  CalFlagFuncArg *largs = (CalFlagFuncArg*)arg;
  ObitUVCalFlagS *me    = largs->in;
  ObitUVDesc *desc      = largs->desc;
  olong      first      = largs->first;
  olong      last       = largs->last;
  olong      iSubA      = largs->iSubA;
  olong      kbase      = largs->kbase;
  olong      SourID     = largs->SourID;
  olong      FQID       = largs->FQID;
  olong      ant1       = largs->ant1;
  olong      ant2       = largs->ant2;
  ofloat     *visIn     = largs->visIn;

  /* local */
  olong iflag, jpoln, jif, jchan, index, limf1, limf2, limc1, limc2;
  olong flga, ipolpt, stadd, incs, incf, incif;

  /* Anything to do? */
  if (last<0) goto finish;
  if (last<first) goto finish;

  /* Increments of things */
  if (desc->inaxes[desc->jlocc]==3) { /* desc correct, complex dim 3 */
    incs  = desc->incs;
    incf  = desc->incf;
    incif = desc->incif;
  } else {  /* multiply by 3 */
    incs  = desc->incs*3;
    incf  = desc->incf*3;
    incif = desc->incif*3;
  }

  /* loop thru flagging criteria */
  ipolpt = abs(me->stoke0)-1;
  if (me->stoke0<-4) ipolpt = abs(me->stoke0)-5;  /* Linear poln */
  for (iflag=first; iflag<=last; iflag++) { /* loop 500 */
 
    /* check antenna */
    flga = me->flagAnt[iflag];
    if ((flga != 0)  &&  (flga != ant1)  &&  (flga != ant2)) continue;

    /* check baseline */
    if ((me->flagBase[iflag] != 0)  &&  (me->flagBase[iflag] != kbase)) continue;

   /* check source */
    if ((me->flagSour[iflag] != SourID)  &&  (me->flagSour[iflag] != 0)  && 
	(SourID != 0)) continue;

    /* check subarray */
    if ((me->flagSubA[iflag] > 0)  &&  (me->flagSubA[iflag] != iSubA)) continue;

    /* check freqid. */
    if ((desc->ilocfq>=0) && (me->flagFQID[iflag] > 0)  &&  
	(me->flagFQID[iflag] != FQID)) continue;

    /* some data to be flagged - set limits */
    limf1 = me->flagBIF[iflag];
    limf2 = me->flagEIF[iflag];
    limc1 = me->flagBChan[iflag];
    limc2 = me->flagEChan[iflag];

    /* loop over polarizations */
    for (jpoln= 0; jpoln<me->numStok; jpoln++) { /* loop 400 */

      if (me->flagPol[iflag*4+jpoln+ipolpt]) { /* Flagged polarization? */
	stadd = jpoln * incs;

	/* loop over IF */
	for (jif= limf1; jif<=limf2; jif++) { /* loop 300 */
	  index = stadd + (jif-1) * incif + (limc1-1) * incf;

	  if (limc1 == limc2) {/* single channel */	    
	    visIn[index+2] = - fabs (visIn[index+2]);

	  } else { /* loop over channel */
	    for (jchan= limc1; jchan<=limc2; jchan++) { /* loop 200 */
	      visIn[index+2] = - fabs (visIn[index+2]);
	      index += incf;
	    } /* end loop over channels L200: */;
	  } 
	} /* end loop  over IF L300: */;
      }
    } /* end loop over Stokes L400: */;
  } /* end loop over flagging entries L500: */;

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadCalFlag */

/**
 * Flag all visibilities in a record
 * \param desc       ObitUVDesc for data
 * \param visIn       visibility as an array of floats
 */
static void CalFlagAll (ObitUVDesc *desc, ofloat *visIn)
{
  olong i;

  for (i=0; i<desc->ncorr; i++) {
    visIn[i*3+2] = - fabs(visIn[i*3+2]);
  }
} /*  end CalFlagAll */

/**
 * Convert a time as time in days to a printable string
 * \param    time  Beginning time, end time in days
 * \msgBuff  Human readable string as "dd/hh:mm:ss.s"
 *           must be allocated at least 13 characters
 */
static void T2String (ofloat time, gchar *msgBuf)
{
  ofloat rtemp, rt1;
  olong   id1, it1, it2;

  id1 = time;
  rtemp = 24.0 * (time - id1);
  it1 = rtemp;
  it1 = MIN (23, it1);
  rtemp = (rtemp - it1)*60.0;
  it2 = rtemp;
  it2 = MIN (59, it2);
  rt1 = (rtemp - it2)*60.0;
  g_snprintf (msgBuf, 30, "%2.2d/%2.2d:%2.2d:%5.2f",
	      id1, it1, it2, rt1);
} /* end of routine T2String   */ 

