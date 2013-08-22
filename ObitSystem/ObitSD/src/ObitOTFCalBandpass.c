/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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

#include "ObitOTFCalBandpassDef.h"
#include "ObitOTFCalBandpass.h"
#include "ObitOTFDesc.h"
#include "ObitOTFSel.h"
#include "ObitTableOTFBP.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFCalBandpass.c
 * ObitOTFCal utilities for applying bandpass calibration to OTF data 
 */

/*-------------------Private function prototypes-------------------*/
/** Private:  Create structure for bandpass Calibration. */
static ObitOTFCalBandpassS* newObitOTFCalBandpassS (ObitOTFCal *in);

/** Private:  Read calibration into the internal arrays. */
static void ObitOTFCalBandpassReadBP (ObitOTFCalBandpassS *in, ObitErr *err);
	
/*----------------------Public functions---------------------------*/
/**
 * Initialize structures for polarization calibration .
 * \param in   Bandpass Object.
 * \param sel  Data selector.
 * \param desc Data descriptor.
 * \param err  ObitError stack.
 */
void ObitOTFCalBandpassInit (ObitOTFCal *in, ObitOTFSel *sel, ObitOTFDesc *desc, 
			     ObitErr *err)
{
  ObitIOCode retCode;
  ObitOTFCalBandpassS *me;
  olong size;
  gchar *routine="ObitOTFCalBandpassInit";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFCalIsA(in));

  in->doBP = sel->doBP;
  if (!in->doBP) return;

  in->bandpassCal = newObitOTFCalBandpassS(in);

  /* pointer to calibration structure */
  me = in->bandpassCal;

  /* Copy Selector information */
  /* me->doBPWt      = sel->doCalWt; Seems better wothout - this is what AIPS does */
  me->doBPWt     = FALSE;
  me->bChan      = sel->startChann;
  me->eChan      = sel->startChann + sel->numberChann - 1;
  me->numFeed    = sel->numberFeed;

  /* Copy Cal information */
  me->BPTable     = ObitTableOTFBPRef(in->BPCalTable);
  me->LastRowRead = 0;

  /* Copy descriptor information */
  me->numChan   = desc->inaxes[desc->jlocf];

  /* Open calibration table  */
  retCode = 
    ObitTableOTFBPOpen ((ObitTableOTFBP*)(me->BPTable), OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);

  /* create row structure */
  me->BPTableRow = (Obit*)newObitTableOTFBPRow((ObitTableOTFBP*)(me->BPTable));
  
  /* table information */
  me->numPol = ((ObitTableOTFBP*)me->BPTable)->numPol;
  me->numRow = ((ObitTableOTFBP*)me->BPTable)->myDesc->nrow;

  /* Allocate calibration arrays */
  /* How big is the calibration table */
  size =(in->numChan) * in->numStok * in->numFeed;
  me->BPApply = g_malloc0(size*sizeof(float));

  /* Polarization offset array PolOff - this allows calibrating 
     both auto and cross pol */
  me->PolOff[0][0] = 0;
  me->PolOff[1][0] = 0;
  me->PolOff[0][1] = me->numChan;
  me->PolOff[1][1] = me->numChan;
  me->PolOff[0][2] = 0;
  me->PolOff[1][2] = me->numChan;
  me->PolOff[0][3] = me->numChan;
  me->PolOff[1][3] = 0;

  /* Read/average BP table */
  ObitOTFCalBandpassReadBP (me, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /*  end ObitOTFCalBandpassInit */

/**
 * Bandpass calibrate data
 * A calibration array entry consists of:
 * Per channel/poln/feed
 * \param in    Bandpass Object.
 * \param time  Time of datum
 * \param DP    Descriptive parameters array.
 * \param recIn 1 record as an array of floats
 * \param err   ObitError stack.
 */
void ObitOTFCalBandpass (ObitOTFCal *in, float time, ofloat *DP, ofloat *recIn, ObitErr *err)
{
  olong   ifeed, ipol, ifreq, ioff, joff, index,nfeed;
  olong   SourID, maxpol, indxa1, indxa2, lenEntryPoln;
  gboolean   sombad, somflg, allflg, smpflg, alpflg, allded;
  gboolean calBad;
  ofloat gwt, g, fblank = ObitMagicF();
  ObitOTFCalBandpassS *me;
  ObitOTFDesc *desc;
  ObitOTFSel *sel;
  /*gchar *routine="ObitOTFCalBandpass";*/

  /* local pointers for structures */
  me   = in->bandpassCal;
  desc = in->myDesc;
  sel  = in->mySel;

  /* Number of feeds */
  nfeed = desc->inaxes[desc->jlocfeed];

  /* Length of entry of all channels with all polarizations */
  lenEntryPoln = MIN (2, MAX (1, in->numStok)) * in->numChan;

  /* Source ID */
  if (desc->iloctar >= 0) SourID = DP[desc->iloctar] + 0.1;
  else SourID = 1;

  /* init. flagged flags */
  allflg = TRUE;
  allded = TRUE;
  alpflg = TRUE;
  somflg = FALSE;
  smpflg = FALSE;

  /* handle 1 solution for both pols - make 0-rel */
  maxpol = MAX (1, MIN(in->numStok, me->numPol))-1;

  /* loop over feed */
  for (ifeed=0; ifeed<nfeed; ifeed++) {
    ioff = ifeed * desc->incfeed;
    
    /* loop over polarization */
    for (ipol=0; ipol<in->numStok; ipol++) {

      /* offset in visibility to this polarization */
      joff = ioff + ipol * desc->incs;

      /* allow auto (ipol=0, 1) and cross (ipol=2,3). */
      indxa1 = ifeed*lenEntryPoln + me->PolOff[0][MIN(ipol,maxpol)];
      indxa2 = ifeed*lenEntryPoln + me->PolOff[1][MIN(ipol,maxpol)];

      /* Initialize calibration */
      g   = 1.0;
      gwt = 1.0;
 
      /* loop over selected channels calibrating. */
      for (ifreq=me->bChan-1; ifreq<me->eChan; ifreq++) {

	/* check if solution valid */
	calBad = (me->BPApply[indxa1+ifreq] == fblank) || 
	         (me->BPApply[indxa2+ifreq] == fblank);

	/* set calibration  */
	if (!calBad) {
	  g = me->BPApply[indxa1+ifreq] * me->BPApply[indxa2+ifreq];
	  
	  /* "weight" calibration */
	  if (me->doBPWt) {
	    gwt = g;
	    if (gwt > 1.0e-10) gwt = 1.0 / gwt;
	  } else {
	    gwt = 1.0;
	  } 
	} else {
	  /* bad calibration - flag data */
	  g   = 0.0;
	  gwt = 0.0;
	}
	
	/* apply calibration */
	index = joff + ifreq * desc->incf;
	/* calibrate data */
	if (recIn[index]!=fblank) {
	  recIn[index]   *= g;
	  recIn[index+1] *= gwt;
	}

	/* keep track of the flagging */
	smpflg = smpflg  ||  (recIn[index+1]  <=  0.0);
	alpflg = alpflg  &&  (recIn[index+1]  <=  0.0);
	somflg = somflg  ||  (gwt  <=  0.0);
	allflg = allflg  &&  (gwt  <=  0.0);
	allded = allded  &&  ((recIn[index+1]  <=  0.0) ||  (gwt  <=  0.0));

      } /* end loop over channels  */
    } /* end loop loop over Stokes  */
  } /* end  loop over Feeds */;

  /* increment counts of the good, bad and the ugly. */
  somflg = somflg  &&  (!allflg);
  smpflg = smpflg  &&  (!alpflg);
  sombad = (somflg || smpflg)  &&  (!allded);
} /* end ObitOTFCalBandpass */

/**
 * Shutdown bandpass calibration.
 * Close any open file and destroy structures.
 * \param in   Bandpass Object.
 * \param err  ObitError stack.
 */
void ObitOTFCalBandpassShutdown (ObitOTFCal *in, ObitErr *err)
{
  ObitOTFCalBandpassS *me;
  ObitIOCode retCode;
  gchar *routine="ObitOTFCalBandpassShutdown";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFCalIsA(in));

  /* Close calibration table  */
  me = in->bandpassCal;
  retCode = ObitTableOTFBPClose ((ObitTableOTFBP*)me->BPTable, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);

  /* delete structure */
  in->bandpassCal = ObitOTFCalBandpassSUnref(in->bandpassCal);
} /*  end ObitOTFCalBandpassShutdown */

/**
 * Destroy structure for bandpass calibration .
 * \param in   Bandpass Object.
 * \return NULL
 */
ObitOTFCalBandpassS*
ObitOTFCalBandpassSUnref (ObitOTFCalBandpassS *in)
{
  if (in==NULL) return in;;
  in->BPTable    = ObitTableOTFBPUnref((ObitTableOTFBP*)in->BPTable);
  in->BPTableRow = ObitTableOTFBPRowUnref((ObitTableOTFBPRow*)in->BPTableRow);
  
  
  if (in->BPApply)      g_free(in->BPApply);  in->BPApply   = NULL;
  /* basic structure */
  g_free (in);

 return NULL;
} /*  end ObitOTFCalBandpassSUnref */

/*---------------Private functions---------------------------*/

/**
 * Create structure for bandpass calibration .
 * \param in   Bandpass Object.
 * \return newly created object.
 */
static ObitOTFCalBandpassS*
newObitOTFCalBandpassS (ObitOTFCal *in)
{
  ObitOTFCalBandpassS* out;

  out = g_malloc0(sizeof(ObitOTFCalBandpassS));

  /* Null pointers */
  out->BPTable      = NULL;
  out->BPTableRow   = NULL;
  out->BPApply      = NULL;

  return out;
} /*  end newObitOTFCalBandpassS */

    
/**
 * Read calibration from BP table.
 * Effectively doBand=1 of interferometry, simple average of all
 * \param in   Bandpass Object.
 * \param err  Error stack for messages and errors.
 */
static void ObitOTFCalBandpassReadBP (ObitOTFCalBandpassS *in, ObitErr *err)
{
  ofloat wt1, fblank = ObitMagicF();
  ofloat *apply=NULL, *sumwt=NULL;
  olong i, istok, ifeed, ichan, indx, lenEntryPoln;
  olong  irow, limit, nchan, nfeed;
  ObitTableOTFBP *BPTable = NULL;
  ObitTableOTFBPRow *BPTableRow = NULL;
  gchar *routine="ObitOTFCalBandpassNewTime";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;

  /* increments, sizes of data elements */
  /* Length of entry of all channels with all polarizations */
  lenEntryPoln = MIN (2, MAX (1, in->numPol)) * in->numChan;
  
  /* Number of feeds */
  nfeed = in->numFeed;

  /* Create/init accumulation arrays */
  apply = g_malloc0(nfeed*lenEntryPoln*sizeof(ofloat));
  sumwt = g_malloc0(nfeed*lenEntryPoln*sizeof(ofloat));

  /* BP table  - set local pointers */
  BPTable = (ObitTableOTFBP*)in->BPTable;
  BPTableRow = (ObitTableOTFBPRow*)in->BPTableRow;

  /* Number of channels in BP table */
  nchan = BPTable->numChan;
 
  /* Read through rows averaging data */
  limit = MAX (1, in->LastRowRead);
  for (i= limit; i<=in->numRow; i++) {
    irow = i;

    /* Read next OTFBP table entry */
    ObitTableOTFBPReadRow (BPTable, irow, BPTableRow, err);
    if (err->error) Obit_traceback_msg (err, routine, "Cal(BP) table");
    if (BPTableRow->status < 0) continue; /* entry flagged? */
    
    /* Averaging all, accumulate in apply and count */
    /* loop over Feeds */
    for (ifeed=0; ifeed<in->numFeed; ifeed++) {
      /* Loop over poln (1 or 2) */
      for (istok=0; istok<in->numPol; istok++) {
	indx = lenEntryPoln * ifeed + istok*nchan;
	wt1  = BPTableRow->wt[indx+ichan];
	/* loop over channels */
	for (ichan=0; ichan<in->numChan; ichan++) {
	  wt1 = BPTableRow->wt[indx+ichan];
	  if ((wt1>0.0) &&  (BPTableRow->mult[indx+ichan]!=fblank)) {
	    apply[indx+ichan] += BPTableRow->mult[indx+ichan] * wt1;
	    sumwt[indx+ichan] += wt1;
	  }
	} /* end channel loop */
      }  /* end stokes loop */
    } /* end feed loop */
  } /* end loop over table */

  /* Normalize to calibration array, 
     take sqrt to allow calibration of X pol data */
  for (i=0; i<nfeed*lenEntryPoln; i++) {
    if ((sumwt[i]>0.0) && (apply[i]!=0.0)) 
      in->BPApply[i] = sqrt (fabs(apply[i] / sumwt[i]));
    else                                
      in->BPApply[i] = fblank;
  }

  /* finished file using all entries */
  in->LastRowRead = in->numRow + 1;

  /* Cleanup */
  if (sumwt) g_free(sumwt);
  if (apply) g_free(apply);
  
} /* end ObitOTFCalBandpassReadBP */
  
