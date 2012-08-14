/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2012                                          */
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

#include <sys/types.h>
#include <time.h>
#include "ObitUVReweight.h"
/*----------------------Private functions---------------------------*/
/** Private: qsort ofloat comparison */
static int compare_gfloat  (const void* arg1,  const void* arg2);


/*----------------------Public functions---------------------------*/
/**
 * Reweight data using the inverse variance of the amplitude in each 
 * data channel.
 * The reweighting is done independently in each time interval defined 
 * by timeAvg. 
 * In each interval, the weights are clipped above the 95 percentile.
 * Need at least 3 samples per interval or the median weight is used.
 * Control parameters on info member of inUV:
 * \li "timeAvg" OBIT_float  (1,1,1) Time interval over which to determine 
 *               reweighting (min) [def = 1 min.]
 *               NB: this should be at least 2 integrations.
 * \li "timeOff" OBIT_float  (1,1,1) Time offset (days) to add to output 
 * \li "subAOff" OBIT_float  (1,1,1) Subarray offset to add to output 
 *
 * \param inUV     Input uv data to reweight. 
 * \param outUV    Output, reweighted data
 *                 May be the same as inUV.
 * \param err      Error stack, returns if  error.
 */
void ObitUVReweightDo (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gboolean doCalSelect;
  olong i, j, k, iwt, jwt, firstVis, startVis, endVis, suba;
  olong countAll;
  olong lastSourceID, curSourceID, lastSubA, lastFQID=-1;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc, *scrDesc;
  ObitUV *scrUV=NULL;
  ofloat timeAvg, lastTime=-1.0, timeOff, subAOff, cbase;
  olong *blLookup=NULL, BIF, BChan;
  olong indx, jndx, nVisPIO, itemp, ant1, ant2;
  olong ncorr, numAnt, numBL, blindx;
  gboolean gotOne, done, incompatible;
  ofloat *acc=NULL, *sortWt=NULL, maxWt, mednWt, *Buffer;
  ofloat startTime, endTime, curTime, rms2, ampl2;
  gchar *routine = "ObitUVReweightDo";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));

  /* Fully instantiate UV data if needed */
  ObitUVFullInstantiate (inUV, TRUE, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  if (outUV) ObitUVFullInstantiate (outUV, TRUE, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);

  /* Get control parameters */
  /* Time interval */
  timeAvg = 1.0;  /* default 1 min */
  ObitInfoListGetTest(inUV->info, "timeAvg", &type, dim, &timeAvg);
  if (timeAvg<=(1.0/60.0)) timeAvg = 1.0;
  timeAvg /= 1440.0;  /* convert to days */
  /* Time offset */
  timeOff = 0.0;
  ObitInfoListGetTest(inUV->info, "timeOff", &type, dim, &timeOff);
  /* Subarray offset */
  subAOff = 0.0;
  ObitInfoListGetTest(inUV->info, "subAOff", &type, dim, &subAOff);

   /* Data Selection */
  BIF = 1;
  ObitInfoListGetTest(inUV->info, "BIF", &type, dim, &BIF);
  if (BIF<1) BIF = 1;
  BChan = 1;
  ObitInfoListGetTest(inUV->info, "BChan", &type, dim, &BChan);
  if (BChan<1) BChan = 1;

 /* Set number of vis per read to 1 */
  nVisPIO = 1;
  ObitInfoListGetTest(inUV->info, "nVisPIO", &type, dim, &nVisPIO);
  itemp = 1;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV->info, "nVisPIO", OBIT_long, dim, &itemp);
  ObitInfoListAlwaysPut(outUV->info, "nVisPIO", OBIT_long, dim, &itemp);

  /* Selection of input? */
  doCalSelect = TRUE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, inUV->name);
  inDesc  = inUV->myDesc;  /* Get descriptor */
  outDesc = outUV->myDesc;

  /* Check compatability between inUV, outUV */
  incompatible = (inDesc->ncorr!=outDesc->ncorr);
  incompatible = incompatible || (inDesc->jlocs!=outDesc->jlocs);
  incompatible = incompatible || (inDesc->jlocf!=outDesc->jlocf);
  incompatible = incompatible || (inDesc->jlocif!=outDesc->jlocif);
  if (incompatible) {
     Obit_log_error(err, OBIT_Error,"%s inUV and outUV have incompatible structures",
		   routine);
      return ;
 }

  /* Check position - same equinox and within 1 mas */
  incompatible = fabs(inDesc->equinox-outDesc->equinox) > 0.01;
  incompatible = incompatible || 
    fabs(inDesc->crval[inDesc->jlocr]-outDesc->crval[inDesc->jlocr]) > 0.001/3600;
  incompatible = incompatible || 
    fabs(inDesc->crval[inDesc->jlocd]-outDesc->crval[inDesc->jlocd]) > 0.001/3600;
  if (incompatible) {
     Obit_log_error(err, OBIT_Error,"%s inUV and outUV have incompatible positions",
		   routine);
      return ;
 }

  /* Check frequency - within 1 Hz */
  incompatible = fabs(inDesc->freq-outDesc->freq) > 1.0;
  if (incompatible) {
     Obit_log_error(err, OBIT_Error,"%s inUV and outUV have incompatible frequencies",
		   routine);
      return ;
 }

  /* Copy input to scratch file to apply editing, selection, flagging */
  iretCode = ObitUVClose (inUV, err);
  scrUV = newObitUVScratch (outUV, err);
  scrUV = ObitUVCopy (inUV, scrUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* Switch in and scr */
  iretCode = ObitUVOpen (scrUV, OBIT_IO_ReadOnly, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, scrUV->name);
  scrDesc = scrUV->myDesc;

  /* Output -
     use same data buffer on input and output 
     so don't assign buffer for output */
  if (outUV->buffer) ObitIOFreeBuffer(outUV->buffer); /* free existing */
  outUV->buffer = NULL;
  outUV->bufferSize = -1;
  oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, outUV->name);
  outDesc->firstVis = outDesc->nvis+1; /* Write to end */

  /* Allocate arrays */
  suba    = 1;
  numAnt  = inUV->myDesc->numAnt[suba-1];/* actually highest antenna number */
  numBL   = (numAnt*(numAnt-1))/2;
  ncorr   = inUV->myDesc->ncorr;
  /* acc index = type + corr * (6) + BL * (3*ncorr)  where BL = 0-rel baseline index */
  acc      = g_malloc0 (ncorr * 3 * numBL * sizeof(ofloat));
  sortWt   = g_malloc0 (ncorr * numBL * sizeof(ofloat));

  /* Baseline tables */
  blLookup = g_malloc0 (numAnt*sizeof(olong));
  blLookup[0] = 0;
  k = 0;
  for (i=1; i<numAnt; i++) {
    blLookup[i] = blLookup[i-1] + numAnt-i;
  }

  /* Initialize things */
  startTime = -1.0e20;
  endTime   =  1.0e20;
  lastSourceID = -1;
  curSourceID  = 0;
  /* Single source selected? */
  if (inUV->mySel->sources) curSourceID = inUV->mySel->sources[0];
  lastSubA     = 0;
  countAll     = 0;
  for (i=0; i<3*ncorr*numBL; i++) acc[i] = 0.0;

  Buffer = scrUV->buffer;

  /* Loop over intervals */
  done   = FALSE;
  gotOne = FALSE;
  startVis  = 1;
  endVis    = firstVis;
  while (!done) {
    
    /* we're in business - loop through data - one vis per read */
    while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
      if ((!gotOne) || (scrDesc->numVisBuff<=0)) { /* need to read new record? */
	iretCode = ObitUVRead (scrUV, NULL, err);
	if (err->error) goto cleanup;
	/*if (iretCode!=OBIT_IO_OK) break;*/
	firstVis = scrDesc->firstVis;
      }
      
      /* Are we there yet??? */
      done = (scrDesc->firstVis >= scrDesc->nvis) || 
	(iretCode==OBIT_IO_EOF);
      if (done && (startTime>0.0)) goto process;

      /* Make sure valid data found */
      if (scrDesc->numVisBuff<=0) continue;
      
      gotOne = FALSE;
      /* Time */
      curTime = Buffer[scrDesc->iloct];
      if (scrDesc->ilocsu>=0) curSourceID = Buffer[scrDesc->ilocsu];
      if (startTime < -1000.0) {  /* Set time window etc. if needed */
	startTime = curTime;
	lastTime  = curTime;
	endTime   = startTime +  timeAvg;
	lastSourceID = curSourceID;
      }

      /* Still in current interval/source? */
      if ((curTime<endTime) && (curSourceID == lastSourceID) && 
	  (scrDesc->firstVis<=scrDesc->nvis) && (iretCode==OBIT_IO_OK)) {
	/* accumulate statistics */
	cbase = Buffer[scrDesc->ilocb]; /* Baseline */
	ant1 = (cbase / 256.0) + 0.001;
	ant2 = (cbase - ant1 * 256) + 0.001;
	lastSubA = (olong)(100.0 * (cbase -  ant1 * 256 - ant2) + 1.5);
	/* Baseline index this assumes a1<a2 always - ignore auto correlations */
	if (ant1!=ant2) {
	  blindx =  blLookup[ant1-1] + ant2-ant1-1;
	  if (scrDesc->ilocfq>=0) lastFQID = (olong)(Buffer[scrDesc->ilocfq]+0.5);
	  else lastFQID = 0;
	  lastTime = curTime;
	  
	  /* Accumulate
	     (1,*) =  count 
	     (2,*) =  sum amp 
	     (3,*) =  sum amp**2 then variance */
	  indx = scrDesc->nrparm; /* offset of start of vis data */
	  for (i=0; i<ncorr; i++) { /* loop 120 */
	    if (Buffer[indx+2] > 0.0) {
	      jndx = i*3 + blindx*3*ncorr;
	      acc[jndx]   += 1.0;
	      ampl2 = Buffer[indx]*Buffer[indx] +  Buffer[indx+1]*Buffer[indx+1];
	      acc[jndx+1] += sqrt(ampl2);
	      acc[jndx+2] += ampl2;
	    } 
	    indx += 3;
	  } /* end loop  L120 over correlations */;
	} /* end only cross correlations */
      } else {  /* process interval */
	
      process:
	/* Now have the next record in the IO Buffer */
	if (iretCode==OBIT_IO_OK) gotOne = TRUE;
	endVis = scrDesc->firstVis;  /* First vis of next block */
	    
	/* Get 1/Variance */
	iwt = 0;
	for (i=0; i<numBL; i++) { /* loop 170 */
	  for (j=0; j<ncorr; j++) { /* loop 160 */
 	    jndx = j*3 + i*3*ncorr;
	    if (acc[jndx] > 2.1) {
	      acc[jndx+1] /= acc[jndx];  /* Average amp */
	      acc[jndx+2] /= acc[jndx];  /* Avg sum amp*amp */
	      /* Variance */
	      rms2 = (acc[jndx+2] - acc[jndx+1]*acc[jndx+1]);
	      rms2 = fabs (rms2);
	      if (rms2>0.0) acc[jndx+2] =  1.0 / rms2;
	      else          acc[jndx+2] =  0.0;
	    } else          acc[jndx+2] = -1.0;  /* Insufficient info */
	    sortWt[iwt++] = acc[jndx+2];  /* List of weights for sorting */
	  } /* end loop  L160: */
	} /* end loop  L170: */

	/* Sort weights */
	qsort ((void*)sortWt, iwt, sizeof(ofloat), compare_gfloat);
	/* Get median or lowest valid */
	mednWt = sortWt[iwt/2];         /* median weight */
	if (mednWt<0.0) {
	  mednWt = 1.0;  /* If all else fails */
	  for (jwt=iwt/2; jwt<iwt; jwt++) {
	    if (sortWt[jwt]>0.0) {mednWt = sortWt[jwt]; break;}
	  }
	}
	/* clip upper 10 % of range */
	maxWt  = sortWt[iwt*95/100];    /* maximum weight */
	if (maxWt<0.0) maxWt = mednWt;

	/* Reinit input to startVis */
	iretCode = ObitUVIOReset (scrUV, startVis, err);
	if (err->error) goto cleanup;

	/* Copy / reweight data */
	while (scrDesc->firstVis<=endVis) {
	  iretCode = ObitUVRead (scrUV, NULL, err);
	  if (err->error) goto cleanup;
	  if (scrDesc->firstVis==endVis) {gotOne=TRUE; break;} /* Finished */
	  if (scrDesc->numVisBuff<=0) continue;  /* any valid data? */

	  /* reweight */
	  countAll++;  /* count */
	  cbase = Buffer[scrDesc->ilocb]; /* Baseline */
	  ant1 = (cbase / 256.0) + 0.001;
	  ant2 = (cbase - ant1 * 256) + 0.001;
	  if (ant1!=ant2) {
	    blindx =  blLookup[ant1-1] + ant2-ant1-1;
	    indx = scrDesc->nrparm; /* offset of start of vis data */
	    for (i=0; i<ncorr; i++) {
	      if (Buffer[indx+2] > 0.0) {
		jndx = i*3 + blindx*3*ncorr;
		if (acc[jndx+2]>=0) Buffer[indx+2] = MIN (maxWt, acc[jndx+2]);
		else                Buffer[indx+2] = mednWt;
	      } 
	      indx += 3;
	    } /* end loop over correlations */;
	  } /* end if cross correlation */

	  /* Modify time, subarray */
	  Buffer[scrDesc->iloct] += timeOff;
	  Buffer[scrDesc->ilocb] += subAOff;

	  /* How many? */
	  outDesc->numVisBuff = scrDesc->numVisBuff;
	  
	  /* Write output */
	  iretCode = ObitUVWrite (outUV, scrUV->buffer, err);
	} /* end loop over data reweighting */

	/* Are we there yet??? */
	done = (scrDesc->firstVis >= scrDesc->nvis) || 
	  (iretCode==OBIT_IO_EOF);

	/* Reinitialize things */
	startTime = -1.0e20;
	endTime   =  1.0e20;
	startVis  = scrDesc->firstVis;
	for (i=0; i<3*ncorr*numBL; i++) acc[i] = 0.0;

      } /* end process interval */
      
    } /* end loop processing data */
    if (done) break;
  } /* end loop over intervals */
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) goto cleanup;

  /* Reset number of vis per read to original value */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV->info, "nVisPIO", OBIT_long, dim, &nVisPIO);
  ObitInfoListAlwaysPut(outUV->info, "nVisPIO", OBIT_long, dim, &nVisPIO);

  /* Inform user */
  Obit_log_error(err, OBIT_InfoErr, "Appended %d reweighted entries", countAll);

  /* Cleanup */
 cleanup:  
  /* close uv files */
  iretCode = ObitUVClose (scrUV, err);
  oretCode = ObitUVClose (outUV, err);
  scrUV = ObitUVUnref(scrUV);  /* Delete scratch file */
  
  /* Deallocate arrays */
  if (acc)      g_free(acc);
  if (sortWt)   g_free(sortWt);
  if (blLookup) g_free(blLookup);
  if (err->error)  Obit_traceback_msg (err, routine, inUV->name);

  return;
} /* end  ObitUVReweightDo */

/*----------------------Private functions---------------------------*/
/**
 * ofloat comparison of two arguments
 * \param arg1 first value to compare
 * \param arg2 second value to compare
 * \return negative if arg1 is less than arg2, zero if equal
 *  and positive if arg1 is greater than arg2.
 */
static int compare_gfloat  (const void* arg1,  const void* arg2)
{
  int out = 0;
  ofloat larg1, larg2;

  larg1 = *(ofloat*)arg1;
  larg2 = *(ofloat*)arg2;
  if (larg1<larg2) out = -1;
  else if (larg1>larg2) out = 1;
  return out;
} /* end compare_gfloat */
