/* $Id$  */
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

#include "ObitOTFCalFlagDef.h"
#include "ObitOTFCalFlag.h"
#include "ObitOTFDesc.h"
#include "ObitOTFSel.h"
#include "ObitTableOTFFlag.h"
#include "ObitTableUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFCalFlag.c
 * ObitOTFCal utilities for applying flagging to uv data 
 */

/*-------------------Private function prototypes-------------------*/
/** Private:  Create structure for flagging. */
static ObitOTFCalFlagS* newObitOTFCalFlagS (ObitOTFCal *in);

/** Private: Update flagging arrays. */
static void ObitOTFCalFlagUpdate (ObitOTFCalFlagS *in, ObitOTFSel *sel, ofloat time,
					ObitErr *err);

/*----------------------Public functions---------------------------*/
/**
 * Initialize structures for Flagging.
 * \param in   Flag Object.
 * \param sel  Data selector.
 * \param desc Data descriptor.
 * \param err  ObitError stack.
 */
void ObitOTFCalFlagInit (ObitOTFCal *in, ObitOTFSel *sel, ObitOTFDesc *desc, 
			 ObitErr *err)
{
  ObitIOCode retCode;
  ObitOTFCalFlagS *me;
  ObitTableOTFFlag *FGTab;
  gchar *routine="ObitOTFCalFlagInit";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFCalIsA(in));

  in->doFlag = sel->doFlag;
  if (!in->doFlag) return;

  in->flag = newObitOTFCalFlagS(in);

  /* pointer to flagging structure */
  me = in->flag;

  /* Copy Selector information */

  /* Copy Cal information */
  me->FlagTable     = ObitTableRef(in->FlagTable);
  me->LastRowRead = 0;
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
  me->numFeed   = desc->inaxes[desc->jlocfeed];
  me->numChan   = desc->inaxes[desc->jlocf];

  /* Make sure sorted in time order */
  FGTab = (ObitTableOTFFlag*)me->FlagTable;
  if (FGTab->myDesc->sort[0]!=(FGTab->TimeRangeOff+1)) {
    retCode = ObitTableUtilSort ((ObitTable*)FGTab, "TIME RANGE ", FALSE, err); 
    if (err->error) Obit_traceback_msg (err, routine, FGTab->name);
  }

  /* Open Flagging table  */
  retCode = 
    ObitTableOTFFlagOpen (FGTab, OBIT_IO_ReadWrite, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);
  
  /* create row structure */
  me->FlagTableRow = (Obit*)newObitTableOTFFlagRow(FGTab);
  
  /* table information */
  me->numRow = FGTab->myDesc->nrow;

  /* Allocate flagging arrays */
  me->maxFlag     = 5000;
  me->numFlag     = 0;
  me->flagTarget  = g_malloc0(me->maxFlag*sizeof(olong));
  me->flagFeed    = g_malloc0(me->maxFlag*sizeof(olong));
  me->flagBChan   = g_malloc0(me->maxFlag*sizeof(olong));
  me->flagEChan   = g_malloc0(me->maxFlag*sizeof(olong));
  me->flagPol     = g_malloc0(4*me->maxFlag*sizeof(gboolean));
  me->flagEndTime = g_malloc0(me->maxFlag*sizeof(ofloat));

  /* Initial time to trigger update of calibration arrays */
  me->flagTime = -1.0e20;
} /*  end ObitOTFCalFlagInit */

/**
 * Flag a record
 * Adapted from AIPS DATFLG.FOR
 * \param in    Flag Object.
 * \param time  Time of datum
 * \param recIn 1 record as an array of floats, 
 * \param err   ObitError stack.
 */
void ObitOTFCalFlag (ObitOTFCal *in, float time, ofloat *recIn, ObitErr *err)
{
  olong   iflag, ifeed, istoke, ichan, TargetID;
  olong limc1, limc2, soff, foff, coff, incdatawt;
  ObitOTFCalFlagS *me;
  ObitOTFDesc *desc;
  ObitOTFSel *sel;
  ofloat *data, fblank = ObitMagicF();
  gboolean doDataWt;
  gchar *routine="ObitOTFCalFlag";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFCalIsA(in));

  /* local pointers for structures */
  me   = in->flag;
  desc = in->myDesc;
  sel  = in->mySel;
  incdatawt = desc->incdatawt; /* increment in data-wt axis */
  doDataWt = incdatawt>1;      /* Have Data-Wt axis? */

  /* see if new time - update cal. */
  if (time > me->flagTime) {
    ObitOTFCalFlagUpdate (me, sel, time, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* If no currently active flags, nothing to do */
  if (me->numFlag <= 0) return;
  
  /* Target ID */
  if (desc->iloctar >= 0) TargetID = recIn[desc->iloctar] + 0.1;
  else TargetID = 0;

  /* loop thru flagging criteria */
  for (iflag= 0; iflag<me->numFlag; iflag++) {
 
    /* check source */
    if ((me->flagTarget[iflag] != TargetID)  &&  (me->flagTarget[iflag] > 0)  && 
	(TargetID != 0)) continue;

    /* some data to be flagged - set limits */
    limc1 = me->flagBChan[iflag];
    limc2 = me->flagEChan[iflag];

    data = &recIn[desc->ilocdata]; /* Data array pointer */

    /* loop over data  */
    for (istoke=0; istoke<me->numStok; istoke++) {
      if (me->flagPol[iflag*4+istoke]) { /* flag this one? */
	soff = istoke * desc->incs; /* offset in data */
	
	/* Loop over feeds */
	for (ifeed=1; ifeed<=me->numFeed; ifeed++) { /* flag this one */
	  if ((ifeed==me->flagFeed[iflag]) || (me->flagFeed[iflag]<=0)) {
	    foff = soff + (ifeed-1) * desc->incfeed;
	    
	    /* Channel loop */
	    for (ichan=limc1; ichan<=limc2; ichan++) {
	      coff = foff + (ichan-1) * desc->incf;
	      data[coff] = fblank;
	      if (doDataWt) data[coff+1] = 0.0;
	      
	    } /* end channel loop */
	  } /* end if feed selected for flagging */
	} /* end Feed loop */
      } /* end if Stokes selected for flagging */
    } /* end Stokes loop */
  } /* end loop over flagging entries */;
} /* end ObitOTFCalFlag */


/**
 * Shutdown baseline dependent Calibrate.
 * Close any open file and destroy structures.
 * \param in   Flag Object.
 * \param err  ObitError stack.
 */
void ObitOTFCalFlagShutdown (ObitOTFCal *in, ObitErr *err)
{
  ObitOTFCalFlagS *me;
  ObitIOCode retCode;
  gchar *routine="ObitOTFCalFlagShutdown";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFCalIsA(in));

  /* Close calibration table  */
  me = in->flag;
  retCode = ObitTableOTFFlagClose ((ObitTableOTFFlag*)me->FlagTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Release row structure  */
  me->FlagTableRow = ObitTableOTFFlagRowUnref(me->FlagTableRow);

  /* delete structure */
  in->flag = ObitOTFCalFlagSUnref(in->flag);
} /*  end ObitOTFCalFlagShutdown */

/**
 * Destroy structure for Flagging .
 * \param in   Flag Object.
 * \return NULL
 */
ObitOTFCalFlagS*
ObitOTFCalFlagSUnref (ObitOTFCalFlagS *in)
{
  if (in==NULL) return in;
  in->FlagTable    = ObitTableOTFFlagUnref((ObitTableOTFFlag*)in->FlagTable);
  in->FlagTableRow = ObitTableOTFFlagRowUnref((ObitTableOTFFlagRow*)in->FlagTableRow);
  if (in->flagTarget)  g_free(in->flagTarget);  in->flagTarget  = NULL;
  if (in->flagFeed)    g_free(in->flagFeed);    in->flagFeed    = NULL;
  if (in->flagBChan)   g_free(in->flagBChan);   in->flagBChan   = NULL;
  if (in->flagEChan)   g_free(in->flagEChan);   in->flagEChan   = NULL;
  if (in->flagPol)     g_free(in->flagPol);     in->flagPol     = NULL;
  if (in->flagEndTime) g_free(in->flagEndTime); in->flagEndTime = NULL;

  /* basic structure */
   g_free (in);

  return NULL;
} /*  end ObitOTFCalFlagSUnref */

/*---------------Private functions---------------------------*/

/**
 * Create structure for baseline dependent Calibrate .
 * \param in   Flag Object.
 * \return newly created object.
 */
static ObitOTFCalFlagS*
newObitOTFCalFlagS (ObitOTFCal *in)
{
  ObitOTFCalFlagS* out;

  out = g_malloc0(sizeof(ObitOTFCalFlagS));

  /* Null pointers */
  out->FlagTable    = NULL;
  out->FlagTableRow = NULL;
  out->flagTarget   = NULL;
  out->flagFeed     = NULL;
  out->flagPol      = NULL;
  out->flagEndTime  = NULL;
  
  return out;
} /*  end newObitOTFCalFlagS */

/**
 * Update ObitOTFCalFlag calibration tables for time time.
 * Adapted from AIPS NXTFLG.FOR
 * \param in   Flag Object.
 * \param sel   OTF selector object
 * \param time desired time in days
 * \param err  Error stack for messages and errors.
 */
static void ObitOTFCalFlagUpdate (ObitOTFCalFlagS *in, ObitOTFSel *sel, ofloat time,
				 ObitErr *err)
{
  olong   ndrop, limit, it;
  olong irow,  i, limit4;
  gboolean done;
  ObitIOCode retCode;
  ObitTableOTFFlag *FlagTable = NULL;
  ObitTableOTFFlagRow *FlagTableRow = NULL;
  gchar *routine="ObitOTFCalFlagUpdate";
 
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert(in!=NULL);
  g_assert(ObitOTFSelIsA(sel));

  /* New time for flags */
  in->flagTime = time;
 
  /* check if any flags expired. */
  done = (in->numFlag <= 0);
  while (!done) {
  
    /* find highest number expired flag */
    ndrop = -1;
    for (i= 0; i<in->numFlag; i++) { /* */
      if (in->flagEndTime[i] < time) ndrop = i;
    } /* end loop */;

    /* see if any to be dropped. */
    done = (ndrop < 0);
    if (done) break;

    /* compress, dropping flag. */
    if ((ndrop+1) < in->numFlag) {  /* This is not the last one the list */
      limit = ndrop;
      for (i= limit; i<in->numFlag-1; i++) {
	in->flagEndTime[i] = in->flagEndTime[i+1];
	in->flagTarget[i]  = in->flagTarget[i+1];
	in->flagFeed[i]    = in->flagFeed[i+1];
	in->flagBChan[i]   = in->flagBChan[i+1];
	in->flagEChan[i]   = in->flagEChan[i+1];
	in->flagPol[i*4]   = in->flagPol[4*(i+1)];
	in->flagPol[i*4+1] = in->flagPol[4*(i+1)+1];
	in->flagPol[i*4+2] = in->flagPol[4*(i+1)+2];
	in->flagPol[i*4+3] = in->flagPol[4*(i+1)+3];
      } /* end loop */;
    } 
    in->numFlag--;
  } /* end loop deleting expired entries */

  /* check if list exhausted */
  if (in->LastRowRead > in->numRow) return;

  /* Flag table - set local pointers */
  FlagTable = (ObitTableOTFFlag*)in->FlagTable;
  FlagTableRow = (ObitTableOTFFlagRow*)in->FlagTableRow;

  /* Find next valid flag. Loop through records */
  limit4 = MAX (1, in->LastRowRead+1);
  for (i= limit4; i<=in->numRow; i++) { 
    irow = i;
    retCode = ObitTableOTFFlagReadRow (FlagTable, irow, FlagTableRow, err);
    if (err->error) Obit_traceback_msg (err, routine, "Flag Table");
    if (FlagTableRow->status < 0) continue; /* entry flagged? */
  
    if (time < FlagTableRow->TimeRange[0]) return;
    if (time > FlagTableRow->TimeRange[1]) continue;
  
    /* Check that starting channel is in range */
    if ((FlagTableRow->chans[0] > 0) && (FlagTableRow->chans[0] > in->numChan)) continue;
    
    /* Is this source wanted */
    if (!ObitOTFSelWantTarget(sel, (olong)FlagTableRow->TargetID)) continue;
    
    /* Must want this one - save it */
    in->LastRowRead = i;
    in->numFlag = in->numFlag + 1;
    /* check if too big */
    if (in->numFlag > in->maxFlag) {
      Obit_log_error(err, OBIT_Error, 
		     "ERROR: Exceeded limit of %d simultaneous flags for %s", 
		     in->maxFlag, sel->name);
      return;
    }
    
    /* fill in tables */
    it = in->numFlag - 1;
    in->flagEndTime[it] = FlagTableRow->TimeRange[1];
    in->flagTarget[it]  = FlagTableRow->TargetID;
    in->flagFeed[it]    = FlagTableRow->Feed;
    in->flagBChan[it] = FlagTableRow->chans[0];
    in->flagEChan[it] = MIN (in->numChan, FlagTableRow->chans[1]);
    if (in->flagBChan[it] <= 0) in->flagBChan[it] = 1;
    if (in->flagEChan[it] <= 0) in->flagEChan[it] = in->numChan;

    /* Ensure that channel selection is in range */
    in->flagEChan[it] = MIN (in->flagEChan[it], in->numChan);

    /* Stokes flags are packed into a bit array */
    in->flagPol[it*4]   = FlagTableRow->pFlags[0] & 0x1;
    in->flagPol[it*4+1] = FlagTableRow->pFlags[0] & 0x2;
    in->flagPol[it*4+2] = FlagTableRow->pFlags[0] & 0x4;
    in->flagPol[it*4+3] = FlagTableRow->pFlags[0] & 0x8;
    
  }  /* end loop over file */

} /* end ObitOTFCalFlagUpdate */
