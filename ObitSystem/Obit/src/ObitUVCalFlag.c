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

#include <glib.h>
#include "ObitUVCalFlag.h"
#include "ObitUVCalFlagDef.h"
#include "ObitUVDesc.h"
#include "ObitUVSel.h"
#include "ObitTableFG.h"
#include "ObitTableUtil.h"

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
  me->maxFlag     = 100000;
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
  me->flagTime = -1.0e20;
} /*  end ObitUVCalFlagInit */

/**
 * Flag a visibility
 * Adapted from AIPS DATFLG.FOR
 * \param in    Flag Object.
 * \param time  Time of datum
 * \param ant1  first antenna number of baseline
 * \param ant2  second antanna of baseline.
 * \param RP    Random parameters array.
 * \param visIn 1 visibility as an array of floats
 * \param err   ObitError stack.
 */
void ObitUVCalFlag (ObitUVCal *in, float time, olong ant1, olong ant2, 
		    ofloat *RP, ofloat *visIn, ObitErr *err)
{
  olong   iflag, kbase, flga,jif, jchan, FQID, SourID, iSubA;
  olong jpoln, limf1, limf2, limc1, limc2, index, stadd, ipolpt;
  ObitUVCalFlagS *me;
  ObitUVDesc *desc;
  ObitUVSel *sel;
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

  /* loop thru flagging criteria */
  for (iflag= 0; iflag<me->numFlag; iflag++) { /* loop 500 */
 
   /* check source */
    if ((me->flagSour[iflag] != SourID)  &&  (me->flagSour[iflag] != 0)  && 
	(SourID != 0)) continue;

    /* check antenna */
    flga = me->flagAnt[iflag];
    if ((flga != 0)  &&  (flga != ant1)  &&  (flga != ant2)) continue;

    /* check baseline */
    if ((me->flagBase[iflag] != 0)  &&  (me->flagBase[iflag] != kbase)) continue;

    /* check subarray */
    if ((me->flagSubA[iflag] > 0)  &&  (me->flagSubA[iflag] != iSubA)) continue;

    /* check freqid. */
    if ((desc->ilocfq>=0) && (me->flagFQID[iflag] > 0)  &&  
	(me->flagFQID[iflag] != FQID)) continue;

    /* some data to be flagged - set limits */
    limf1 = MAX (1, me->flagBIF[iflag]);
    limf2 = MIN (me->numIF, me->flagEIF[iflag]);
    limc1 = MAX (1, me->flagBChan[iflag]);
    limc2 = MIN (me->numChan,  me->flagEChan[iflag]);

    /* loop over polarizations */
    ipolpt = abs(me->stoke0)-1;
    for (jpoln= 0; jpoln<me->numStok; jpoln++) { /* loop 400 */

      if (me->flagPol[iflag*4+jpoln+ipolpt]) { /* Flagged polarization? */
	stadd = jpoln * desc->incs;

	/* loop over IF */
	for (jif= limf1; jif<=limf2; jif++) { /* loop 300 */
	  index = stadd + (jif-1) * desc->incif + (limc1-1) * desc->incf;

	  if (limc1 == limc2) {/* single channel */	    
	    visIn[index+2] = - fabs (visIn[index+2]);

	  } else { /* loop over channel */
	    for (jchan= limc1; jchan<=limc2; jchan++) { /* loop 200 */
	      visIn[index+2] = - fabs (visIn[index+2]);
	      index += desc->incf;
	    } /* end loop over channels L200: */;
	  } 
	} /* end loop  over IF L300: */;
      }
    } /* end loop over Stokes L400: */;
  } /* end loop over flagging entries L500: */;
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
  olong   ndrop, limit, a1, a2, it;
  olong irow,  i, limit4;
  gboolean done;
  ObitIOCode retCode;
  ObitTableFG *FGTable = NULL;
  ObitTableFGRow *FGTableRow = NULL;
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
    for (i= 0; i<in->numFlag; i++) { /* loop 100 */
      if (in->flagEndTime[i] < time) ndrop = i;
    } /* end loop  L100: */;

    /* see if any to be dropped. */
    done = (ndrop < 0);
    if (done) break;

    /* compress, dropping flag. */
    if ((ndrop+1) < in->numFlag) {  /* This is not the last one the list */
      limit = ndrop;
      for (i= limit; i<in->numFlag-1; i++) { /* loop 150 */
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
      } /* end loop  L150: */;
    } 
    in->numFlag--;
  } /* end loop deleting expired entries */

  /* check if list exhausted */
  if (in->LastRowRead >= in->numRow) return;

  /* FG table - set local pointers */
  FGTable = (ObitTableFG*)in->FGTable;
  FGTableRow = (ObitTableFGRow*)in->FGTableRow;

  /* Find next valid flag. Loop through records */
  limit4 = MAX (1, in->LastRowRead+1);
  for (i= limit4; i<=in->numRow; i++) { /* loop 360 */
    irow = i;
    retCode = ObitTableFGReadRow (FGTable, irow, FGTableRow, err);
    if (err->error) Obit_traceback_msg (err, routine, "FG Table");
    if (FGTableRow->status < 0) continue; /* entry flagged? */
    in->LastRowRead = i-1;  /* may try this one again */
  
    if (time < FGTableRow->TimeRange[0]) return;
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
    in->numFlag = in->numFlag + 1;
    /* check if too big */
    if (in->numFlag > in->maxFlag) {
      Obit_log_error(err, OBIT_Error, 
		     "ERROR: Exceeded limit of %d simultaneous flags for %s", 
		     in->maxFlag, sel->name);
      return;
    }
    
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
      in->flagAnt[it] = FGTableRow->ants[0];
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

    /* Ensure that if and channel selection are in range */
    in->flagEIF[it] = MIN (in->flagEIF[it], in->numIF);
    in->flagEChan[it] = MIN (in->flagEChan[it], in->numChan);

    /* Stokes flags are packed into a bit array */
    in->flagPol[it*4]   = FGTableRow->pFlags[0] & 0x1;
    in->flagPol[it*4+1] = FGTableRow->pFlags[0] & 0x2;
    in->flagPol[it*4+2] = FGTableRow->pFlags[0] & 0x4;
    in->flagPol[it*4+3] = FGTableRow->pFlags[0] & 0x8;
    
  }  /* end loop over file  L360: */;

} /* end ObitUVCalFlagUpdate */
