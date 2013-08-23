/* $id:  $  */
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

#include <time.h>
#include "ObitUtil.h"
#include "ObitOTF.h"
#include "ObitTableOTFIndex.h"
#include "ObitTableOTFTargetUtil.h"
#include "ObitTableOTFFlag.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_multifit_nlin.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFCalUtil.c
 * ObitOTF class data flagging utility function definitions.
 */

/*---------------Private structures----------------*/
/* FD Median calculations threaded function argument */
typedef struct {
  /* ObitThread with restart queue */
  ObitThread *thread;
  /* thread number, <0 -> no threading   */
  olong        ithread;
  /* First (0-rel) Spectrum to process this thread */
  olong        first;
  /* Highest (0-rel) Spectrum to process this thread  */
  olong        last;
  /* Number of Channels */
  olong numChan; 
  /*  Count of entries in each cell of avg, RMS */
  ollong *count;
  /* Average (freq, poln, feed) */
  ofloat *avg;
  /* Amp^2 in RMS out (freq, poln, feed */
  ofloat *RMS;
  /* [out] median alpha sigma for values in avg  */
  ofloat *sigma;
  /* Max. residual flux from running median in sigma allowed */
  ofloat maxRes;
  /* Flag all channels having RMS values > maxRMS[0] of the 
     channel median sigma.[default 6.0] plus maxRMS[1] (default 0.1) 
     of the channel average in quadrature
     if maxRMS[0] > 100.0 no RMS flagging */
  ofloat *maxRMS;
   /* Min. fraction of good times */
  ofloat minGood;
  /* the width of the median window in channels */
  olong widMW;  
  /* Error stack, returns if error. */
  ObitErr* err;   
  /* Work arrays at least the size of numChan */
  ofloat *work, *work2;
} OTFMednFuncArg;

/*---------------Private function prototypes----------------*/
/** Private: Process for FD */
static void EditFDProcess(olong nThread, OTFMednFuncArg** args, ObitErr *err);

/** Private: Threaded FD median/linear editing */
static gpointer ThreadEditFDProcess(gpointer arg);

/** Private: Determine sigma for Median */
static ofloat FDMedianSigma (olong n, ofloat *value, ofloat mean, ofloat alpha);

/*----------------------Public functions---------------------------*/

/**
 * Frequency domain editing of OTF data.  Possibly threaded.
 * Editing is done independently for each data stream.  
 * An average and RMS is determined for each channel  in each timeAvg 
 * period and a spectral baseline is established for the average values
 * using a running median window.  Channels with excessive RMSes or 
 * residual amplitudes are flagged.  Channels with too few samples may 
 * also be edited.
 * Flagging is done by entering the offending data in OTFFlag table flagTab
 * on outOTF.
 * \param inOTF     Name of input uvdata object. 
 *                 Any prior selection and editing is applied.
 * Control parameters on info member of inOTF:
 * \li "flagTab" OBIT_int    (1,1,1) FG table version number [ def. 1]
 * \li "timeAvg" OBIT_float  (1,1,1) Time interval over which to average 
 *               data to be flagged (min) [def = 1 min.]
 *               NB: this should be at least 2 integrations.
 * \li "FDwidMW"  OBIT_int (1,1,1) The width of the median window in channels. 
 *               An odd number (5) is recommended,  default or 0 -> linear baseline
 * \li "FDmaxRMS" OBIT_float (2,1,1) Flag all channels having RMS 
 *               values > maxRMS[0] of the channel median sigma.[default => 6.]
 *               plus maxRMS[1] (default 0.1) of the channel average in quadrature
 * \li "FDmaxRes" OBIT_float (1,1,1) Max. residual flux in sigma.  default => 6.
 * \li "minGood" OBIT_float (1,1,1) Minimum fraction of time samples 
 *               at and below which a channel/interval will be flagged.  
 *               [default 0.25, -1 ->no flagging by fraction of samples.]
 * \param outOTF   OTF data onto which the FG table is to be attached.
 *                 May be the same as inOTF.
 * \param err      Error stack, returns if not empty.
 */
void ObitOTFFlagEditFD (ObitOTF* inOTF, ObitOTF* outOTF, ObitErr* err) 
{
  ObitIOCode iretCode, oretCode;
  ObitTableOTFFlag *tmpFlag=NULL;
  ObitTableOTFFlag *outFlag=NULL;
  ObitTableOTFFlagRow *row=NULL;
  gboolean doCalSelect;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  ObitOTFDesc *inDesc;
  OTFMednFuncArg** args=NULL;
  ObitThread *myThread=NULL;
  olong flagTab, iFGRow, nThread, nSppTh;
  ofloat timeAvg, maxRes, maxRMS[2], minGood;
  olong  numChan, numStok, numFeed, widMW, incf, incs, incfeed;
  olong jndx, indx, kndx, jj, ifeed, istok, ichan;
  olong i, kk, itemp, lastSourceID, curSourceID, irow;
  ollong  countAll, countBad, *count=NULL;
  olong nRecPIO, ver, firstRec, startRec;
  ofloat startTime, endTime, curTime, *Buffer;
  ofloat lastTime=-1.0;
  gboolean done, gotOne;
  /* Accumulators per spectral channel/poln/feed */
  ofloat *sumA=NULL, *sumA2=NULL, *sigma=NULL;
  ofloat sec;
  gchar *tname, reason[25];
  struct tm *lp;
  time_t clock;
  gchar *routine = "ObitOTFFlagEditFD";

  /* error checks */
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));
  g_assert (ObitOTFIsA(outOTF));

  /* Fully instantiate OTF data if needed */
  ObitOTFFullInstantiate (inOTF, TRUE, err);
  if (outOTF) ObitOTFFullInstantiate (outOTF, TRUE, err);
  if (err->error) Obit_traceback_msg (err, routine, outOTF->name);

  /* Get control parameters */
  flagTab = 1;
  ObitInfoListGetTest(inOTF->info, "flagTab", &type, dim, &flagTab);
  /* Time interval */
  timeAvg = 1.0;  /* default 1 min */
  ObitInfoListGetTest(inOTF->info, "timeAvg", &type, dim, &timeAvg);
  if (timeAvg<=(1.0/60.0)) timeAvg = 1.0;
  timeAvg /= 1440.0;  /* convert to days */
  /* Median window size */
  widMW = 31;  
  ObitInfoListGetTest(inOTF->info, "FDwidMW", &type, dim, &widMW);
  /* Maximum average residual */
  maxRes = 6.0;  
  ObitInfoListGetTest(inOTF->info, "FDmaxRes", &type, dim, &maxRes);
  if (maxRes<=0.0) maxRes = 6.0; 
  /* Maximum RMS */
  maxRMS[0] = 6.0; maxRMS[1] = 0.1;
  ObitInfoListGetTest(inOTF->info, "FDmaxRMS", &type, dim, maxRMS);
  if (maxRMS[0]<=0.0) maxRMS[0] = 6.0; 
   /* Min. fraction of good times */
  minGood = 0.25;           /* default 0.25 */
  ObitInfoListGetTest(inOTF->info, "minGood", &type, dim,  &minGood);  
 
 /* Set number of rec per read to 1 */
  nRecPIO = 1;
  ObitInfoListGetTest(inOTF->info, "nRecPIO", &type, dim, &nRecPIO);
  itemp = 1;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inOTF->info, "nRecPIO", OBIT_long, dim, &itemp);

  /* Selection of input? */
  doCalSelect = TRUE;
  ObitInfoListGetTest(inOTF->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* Open input */
  iretCode = ObitOTFOpen (inOTF, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, inOTF->name);
  inDesc  = inOTF->myDesc;  /* Get descriptor */

  /* Create temporary output FG table */
  tname = g_strconcat ("Temp Flag Table for: ", outOTF->name, NULL);
  ver = 0;
  tmpFlag = newObitTableOTFFlagValue(tname, (ObitData*)outOTF, &ver, OBIT_IO_WriteOnly, 
				     err);
  g_free (tname);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
 
  /* Allocate arrays */
  if (inOTF->myDesc->jlocs>=0)    numStok = inOTF->myDesc->inaxes[inOTF->myDesc->jlocs];
  else numStok = 1;
  if (inOTF->myDesc->jlocfeed>=0) numFeed = inOTF->myDesc->inaxes[inOTF->myDesc->jlocfeed];
  else numFeed = 1;
  if (inOTF->myDesc->jlocf>=0)    numChan = inOTF->myDesc->inaxes[inOTF->myDesc->jlocf];
  else numChan = 1;
  incf    = inOTF->myDesc->incf;
  incs    = inOTF->myDesc->incs;
  incfeed = inOTF->myDesc->incfeed;
  /* index in order channel/Stokes/feed */
  sumA     = g_malloc0 (numChan*numStok*numFeed*sizeof(ofloat));
  sumA2    = g_malloc0 (numChan*numStok*numFeed*sizeof(ofloat));
  sigma    = g_malloc0 (numChan*numStok*numFeed*sizeof(ofloat));
  count    = g_malloc0 (numChan*numStok*numFeed*sizeof(ollong));

  /* How many Threads? */
  myThread = newObitThread();
  nThread = MAX (1, ObitThreadNumProc(myThread));
  nThread = MIN (nThread, numStok*numFeed);
  /* Number of spectra per thread */
  nSppTh  = MAX (1, (olong)(0.999+(ofloat)(numFeed*numStok)/nThread));  /* Spectra per thread */
  nSppTh  = MIN (nSppTh, numFeed*numStok);

  /* Create Thread object arrays */
  args = g_malloc0(nThread*sizeof(OTFMednFuncArg*));
  for (i=0; i<nThread; i++) {
    args[i] = g_malloc0(sizeof(OTFMednFuncArg));
    args[i]->thread = myThread;
    if (nThread>1) args[i]->ithread = i;
    else           args[i]->ithread = -1;
    args[i]->first    = i*nSppTh;
    args[i]->last     = MIN(numFeed*numStok-1, args[i]->first + nSppTh-1);
    if (i==(nThread-1)) args[i]->last = numFeed*numStok-1;  /* Make sure all done */
    args[i]->numChan  = numChan;
    args[i]->count    = count;
    args[i]->avg      = sumA;
    args[i]->RMS      = sumA2;
    args[i]->sigma    = sigma;
    args[i]->minGood  = minGood;
    args[i]->maxRMS   = maxRMS;
    args[i]->maxRes   = maxRes;
    args[i]->widMW    = MIN (widMW, numChan-1);
    args[i]->err      = err;
    args[i]->work     = g_malloc(numChan*sizeof(ofloat));   /* Deallocate when done */
    args[i]->work2    = g_malloc(numChan*sizeof(ofloat));   /* Deallocate when done */
  }

  /* Initialize things */
  startTime = -1.0e20;
  endTime   =  1.0e20;
  lastSourceID = -1;
  curSourceID  = 0;
  /* Single source selected? */
  if (inOTF->mySel->targets) curSourceID = inOTF->mySel->targets[0];
  countAll     = 0;
  countBad     = 0;

  Buffer = inOTF->buffer;

  /* Open output table */
  oretCode = ObitTableOTFFlagOpen (tmpFlag, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  
  /* Create Row */
  row = newObitTableOTFFlagRow (tmpFlag);
  
  /* Attach row to output buffer */
  ObitTableOTFFlagSetRow (tmpFlag, row, err);
  if (err->error) goto cleanup;
  
  /* If there are entries in the table, mark it unsorted */
  if (tmpFlag->myDesc->nrow>0) 
    {tmpFlag->myDesc->sort[0]=0; tmpFlag->myDesc->sort[1]=0;}
  
  /* Initialize solution row */
  row->TargetID= 0; 
  row->Feed  = 0; 
  row->TimeRange[0] = -1.0e20; 
  row->TimeRange[1] =  1.0e20; 
  row->chans[0]  = 0; 
  row->chans[1]  = 0; 
  row->pFlags[0] = 0; 
  row->pFlags[1] = 0; 
  row->pFlags[2] = 0; 
  row->pFlags[3] = 0; 
  /* Reason includes time/date */
  /* Get time since 00:00:00 GMT, Jan. 1, 1970 in seconds. */
  time (&clock);
  /* Convert to  broken-down time. */
  lp = localtime (&clock);
  lp->tm_mon++; /* For some bizzare reason, month is 0-rel */
  if (lp->tm_year<1000)  lp->tm_year += 1900; /* full year */
  sec = (ofloat)lp->tm_sec;
  g_snprintf (reason, 25, "EditFD %d/%d/%d %d:%d:%f", 
	      lp->tm_year, lp->tm_mon, lp->tm_mday, 
	      lp->tm_hour, lp->tm_min, sec);
  row->reason    = reason; /* Unique string */
  
  /* Loop over intervals */
  done   = FALSE;
  gotOne = FALSE;
  while (!done) {
    
    /* we're in business - loop through data - one rec per read */
    while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
      if ((!gotOne) || (inOTF->myDesc->numRecBuff<=0)) { /* need to read new record? */
	if (doCalSelect) iretCode = ObitOTFReadSelect (inOTF, inOTF->buffer, err);
	else iretCode = ObitOTFRead (inOTF, inOTF->buffer, err);
	if (err->error) goto cleanup;
	/*if (iretCode!=OBIT_IO_OK) break;*/
	firstRec = inDesc->firstRec;
      }
      
      /* Are we there yet??? */
      done = (inDesc->firstRec >= inDesc->nrecord) || 
	(iretCode==OBIT_IO_EOF);
      if (done && (startTime>0.0)) goto process;

      /* Make sure valid data found */
      if (inOTF->myDesc->numRecBuff<=0) continue;
      
      gotOne = FALSE;
      /* Time */
      curTime = Buffer[inDesc->iloct];
      if (inDesc->iloctar>=0) curSourceID = Buffer[inDesc->iloctar];
      if (startTime < -1000.0) {  /* Set time window etc. if needed */
	startTime = curTime;
	lastTime  = curTime;
	endTime   = startTime +  timeAvg;
	startRec  = firstRec;
	lastSourceID = curSourceID;
      }

      /* Still in current interval/source? */
      if ((curTime<endTime) && (curSourceID == lastSourceID) && 
	  (inDesc->firstRec<=inDesc->nrecord) && (iretCode==OBIT_IO_OK)) {
	/* accumulate statistics */
	lastTime = curTime;
	  
	/* Accumulate
	   count =  count per channel, poln, feed
	   sumA  = sum of amplitudes for channel, poln, feed, then residual
	   sumA2 = sum of amplitudes**2 for channel, poln, feed, then RMS */
	/* Loop over feed */
	for (ifeed=0; ifeed<numFeed; ifeed++) {
	  /* Loop over Stokes */
	  for (istok=0; istok<numStok; istok++) {
	    /* Loop over Channel */
	    for (ichan=0; ichan<numChan; ichan++) {
	      indx = inDesc->ilocdata + ichan*incf + istok*incs + ifeed*incfeed;
	      jndx = ifeed*numChan*numStok + istok*numChan + ichan;
	      if (Buffer[indx+1]>0.0) {
		count[jndx]++;
		sumA[jndx]  += Buffer[indx];
		sumA2[jndx] += Buffer[indx]*Buffer[indx];
	      }
	    } /* end channel loop */
	  } /* end Stokes loop */
	} /* end feed loop */
      } else {  /* process interval */
	
      process:
	/* Now have the next record in the IO Buffer */
	if (iretCode==OBIT_IO_OK) gotOne = TRUE;

	/* Process  */
	EditFDProcess(nThread, args, err); 
	if (err->error) goto cleanup;
	    
	/* Init Flagging table entry */
	row->TargetID     = lastSourceID; 
	row->TimeRange[0] = startTime;
	row->TimeRange[1] = lastTime;
	
	/* Loop over feed */
	for (ifeed=0; ifeed<numFeed; ifeed++) {
	  /* Loop over Stokes */
	  for (istok=0; istok<numStok; istok++) {
	    /* Loop over Channel */
	    for (ichan=0; ichan<numChan; ichan++) {
	      jndx = ifeed*numChan*numStok + istok*numChan + ichan;
	      countAll++;  /* How many examined */
	      if (sumA[jndx]<-9900.0) {
		countBad++;  /* Count flagged interval/correlator */
		
		/* Check for higher number spectral channels in a contigious 
		   block of bad channels and include in same flag */
		  kndx = jndx+1;
		  jj = ichan;
		  for (kk = ichan+1; kk<numChan; kk++) {
		    /* Higher channel number and to be flagged? */
		    if (sumA[kndx]<-9900.0) { /* This one flagged? */
		      jj = kk;          /* This correlator to be included in flag */
		      countBad++;       /* Increment bad count */
		      sumA[kndx] = 0.0; /* don't consider again */
		    } else break;       /* stop looking at higher channel numbers */
		    kndx++;
		  } /* end loop searching for higher bad channels */
		  /* Flag channels ichan+1, jj+1 (1-rel) */
		  row->Feed     = ifeed+1;
		  row->chans[0] = ichan+1; row->chans[1] = jj+1; 
		  if (istok==0) row->pFlags[0] = 1;
		  else if (istok==1) row->pFlags[0] = 2;
		  else if (istok==2) row->pFlags[0] = 4;
		  else if (istok==3) row->pFlags[0] = 8;
		  else row->pFlags[0] = 0;
		  /* Write */
		  iFGRow = -1;
		  oretCode = ObitTableOTFFlagWriteRow (tmpFlag, iFGRow, row, err);
		  if (err->error) goto cleanup;
	      } /* end if flagged */
	    } /* end channel loop */
	  } /* end Stokes loop */
	} /* end feed loop */

	/* Are we there yet??? */
	done = (inDesc->firstRec >= inDesc->nrecord) || 
	  (iretCode==OBIT_IO_EOF);

	/* Reinitialize things */
	startTime = -1.0e20;
	endTime   =  1.0e20;
	for (i=0; i<numChan*numStok*numFeed; i++) {
	  sumA[i]  = 0.0;
	  sumA2[i] = 0.0;
	  count[i] = 0;
	}
      } /* end process interval */
      
    } /* end loop processing data */
    if (done) break;
  } /* end loop over intervals */
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) goto cleanup;
  
  /* close uv file */
 cleanup:
  iretCode = ObitOTFClose (inOTF, err);
  
  /* Copy temp OTFFlag table tmpFlag contents to version flagTab 
     if tmpFlag is version > flagTab */
  if (ver>flagTab) {
    tname = g_strconcat ("Output Flag Table for: ", outOTF->name, NULL);
    ver     = flagTab;
    outFlag = newObitTableOTFFlagValue(tname, (ObitData*)outOTF, &ver, OBIT_IO_ReadWrite, 
				       err);
    g_free (tname);
    
    /* Open output table */
    oretCode = ObitTableOTFFlagOpen (outFlag, OBIT_IO_ReadWrite, err);
    if (err->error) goto cleanup2;
    
    /* Attach row to output buffer */
    ObitTableOTFFlagSetRow (outFlag, row, err);
    if (err->error) goto cleanup2;
    
    /* If there are entries in the table, mark it unsorted */
    if (outFlag->myDesc->nrow>0) 
      {outFlag->myDesc->sort[0]=0; outFlag->myDesc->sort[1]=0;}
    
    /* Append */
    for (irow=1; irow<=tmpFlag->myDesc->nrow; irow++) {
      iretCode = ObitTableOTFFlagReadRow (tmpFlag, irow, row, err);
      iFGRow = -1;
      oretCode = ObitTableOTFFlagWriteRow (outFlag, iFGRow, row, err);
      if (err->error) goto cleanup2;
    } /* End copy */
    
    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "EditFD: Copied %d flags to table %d",
		   tmpFlag->myDesc->nrow, ver);
    /* Close flag tables */
    oretCode = ObitTableOTFFlagClose (outFlag, err);
    oretCode = ObitTableOTFFlagClose (tmpFlag, err);
 
    /* Delete tmpFlag */
    tmpFlag = (ObitTableOTFFlag*)ObitTableZap((ObitTable*)tmpFlag, err);

  } else { /* end copy */
    /* Only one - keep it */
    oretCode = ObitTableOTFFlagClose (tmpFlag, err);
  } /* end keep OTFFlag table */

  /* Reset number of rec per read to original value */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inOTF->info, "nRecPIO", OBIT_long, dim, &nRecPIO);

  /* Cleanup */
  /* Deallocate arrays */
 cleanup2:
  row     = ObitTableOTFFlagRowUnref(row);
  outFlag = ObitTableOTFFlagUnref(outFlag);
  if (sumA)     g_free(sumA);
  if (sumA2)    g_free(sumA2);
  if (count)    g_free(count);
  if (sigma)    g_free(sigma);
  /* Shut down any threading */
  ObitThreadPoolFree (myThread);
  /* Delete Thread object arrays */
  if (args) {
    for (i=0; i<nThread; i++) {
      if (args[i]) {
	if (args[i]->work)  g_free(args[i]->work);
	if (args[i]->work2) g_free(args[i]->work2);
	g_free(args[i]); args[i] = NULL;
      }
    }
    g_free(args); args = NULL;
  }
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Give report */
  Obit_log_error(err, OBIT_InfoErr, "EditFD: flag %8.1lf of %8.1lf rec/interval= %5.1lf percent",
		 (odouble)countBad, (odouble)countAll, 
		 100.0*(odouble)countBad/((odouble)countAll));

  return;
} /* end  ObitOTFEditFD */


/*----------------------Private functions---------------------------*/
/**
 * Do processing for FD editing, possible with threading
 * \param nThread  Number of parallel threads
 * \param args     arguments per thread
 * \param err      Error stack, returns if  error.
 */
static void EditFDProcess(olong nThread, OTFMednFuncArg** args, ObitErr *err)
{
  gboolean OK;
  ObitThreadFunc func=(ObitThreadFunc)ThreadEditFDProcess;
  gchar *routine = "EditFDprocess";

  /* error checks */
  if (err->error) return;

  /* Do operation onarrays possibly with threads */
  OK = ObitThreadIterator (args[0]->thread, nThread, func, (gpointer**)args);

  /* Check for problems */
  if (!OK) {
    Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
  }

} /* end EditFDProcess */

/**
 * Threaded Processing for Frequency domain editing
 * Flag channels with less than or equal minGood of time samples valid
 * \param arg     pointer to structure with:
 * \li first    First (0-rel) Spectrum to process this thread
 * \li last     Highest (0-rel) Spectrum to process this thread
 * \li numChan  Number of spectral channels in count, sumA, RMS
 * \li count    Count of entries in each cell of avg, RMS
 * \li avg      Sum of amplitudes in time interval for freq, IF, poln, baseline
 *                  <=0.0 => no data or flagged
 * \li RMS      Sum of amplitude**2 in time interval for freq, IF, poln, baseline
 * \li sigma    [out] median alpha sigma for values in avg (MW only)
 *                corresponding Stokes RR for a Stokes LL measurement, used with maxV
 * \li widMW    If > 0 the width of the median window in channels.
 *                 If <= 0 use linear baseline and chanMask
 * \li maxRMS   Flag all channels having RMS values > maxRMS[0] of the 
 *              channel median sigma.[default 6.0] plus maxRMS[1] (default 0.1) 
 *              of the channel average in quadrature
 *              if maxRMS[0] > 100.0 no RMS flagging
 * \li maxRes   Max. residual flux in sigma allowed for channels outside 
 *              the baseline fitting regions.
 * \li err      Error stack, returns if not empty.
 */
static gpointer ThreadEditFDProcess(gpointer arg)
{
  /* Get arguments from structure */
  OTFMednFuncArg *largs = (OTFMednFuncArg*)arg;
  olong loSp         = largs->first;
  olong hiSp         = largs->last;
  ofloat *sumA       = largs->avg;
  ofloat *avg        = largs->avg;
  ofloat *sumA2      = largs->RMS;
  ofloat *RMS        = largs->RMS;
  ofloat *sigma      = largs->sigma;
  ollong *count      = largs->count;
  olong  numChan     = largs->numChan; 
  olong  widMW       = largs->widMW;
  ofloat *maxRMS     = largs->maxRMS;
  ofloat maxRes      = largs->maxRes;
  ofloat minGood     = largs->minGood;

  olong jf, jsp, indx, jndx, jj, jjj, cnt, half;
  ollong maxCount, flagCount;
  ofloat aaa, *temp, *temp2;
  odouble sum1, sum2, meanRMS, residualRMS;

  /* error checks */
  if (largs->err->error) goto done;

  /* Loop over spectra averaging */
  for (jsp=loSp; jsp<=hiSp; jsp++) {
    indx = jsp*numChan;
    for (jf=0; jf<numChan; jf++) {
      if (count[indx]>0) {  /* Any data? */
	sumA[indx] /= count[indx];
	if (count[indx]>1) { /* Enough data for RMS? */
	  sumA2[indx] = sqrt((sumA2[indx]/count[indx] - sumA[indx]*sumA[indx])*
			     ((odouble)count[indx])/(count[indx]-1.0));
	} else sumA2[indx] = 0.0;
      } /* end if data */
      indx++;
    } /* end loop over Channel */
  } /* end loop over spectrum */

  /* Subtract running medians from averaged data */
  temp  = largs->work;
  temp2 = largs->work2;
  half  = widMW/2;  /* Half width of window */
  /* Loop over spectra */
  for (jsp=loSp; jsp<=hiSp; jsp++) {
    indx = jsp*numChan;
    for (jf=0; jf<numChan; jf++) temp2[jf] =  avg[indx+jf];
    for (jf=0; jf<numChan; jf++) {
      /* Median window are we at the beginning, end or middle? */
      cnt = 0;
      if (jf<half) { /* Beginning */
	for (jj=0; jj<widMW; jj++) {
	  jjj = jj - jf;  /* If we're no longer at the beginning */
	  if ((count[indx+jjj]>0) && (temp2[jj]>-9900.0)) {  /* Any data? */
	    temp[cnt++] = temp2[jj];
	  }
	}
      } else if (jf>=(numChan-half-1)) { /* End */
	for (jj=numChan-widMW; jj<numChan; jj++) {
	  jjj = jj - jf;  /* If we're no longer at the beginning */
	  if ((count[indx+jjj]>0) && (temp2[jj]>-9900.0)) {  /* Any data? */
	    temp[cnt++] =  temp2[jj];
	  }
	}
      } else { /* Middle */
	for (jj=-half; jj<=half; jj++) {
	  if ((count[indx+jj]>0) && (temp2[jf+jj]>-9900.0)) {  /* Any data? */
	    temp[cnt++] = temp2[jf+jj];
	  }
	}
      }
      /* Get Median and subtract */
      if ((count[indx]>0) && (avg[indx]>-9900.0)) {
	aaa = medianValue (temp, 1, cnt);
	avg[indx] -= aaa;
	sigma[indx] = FDMedianSigma (cnt, temp, aaa, 0.2);
      } else {  /* Datum bad */
	if (avg[indx]>-9900.0) avg[indx]   = 0.0;
	sigma[indx] = 100.0;
      }
      indx++;
    } /* end loop over Channel */
  } /* end loop over spectrum */
  

  /* Decide which to edit */
  temp = largs->work; /* Work array for median */
 
  /* Loop over spectra*/
  for (jsp=loSp; jsp<=hiSp; jsp++) {
    indx = jsp*numChan;
    jndx = 0;
    /* Determine average RMS for spectrum */
    if (maxRMS[0]<100.0) {
      cnt = 0;
      for (jf=0; jf<numChan; jf++) {
	if ((avg[indx+jf]>-9900.0) && (RMS[indx+jf]>0.0)) {
	  temp[cnt++] = RMS[indx+jf];
	}
      }
      
      /* Median rms */
      if (cnt>0) meanRMS = medianValue(temp, 1, cnt);
      else meanRMS = 0.0;
      
      /* RMS flagging */
      for (jf=0; jf<numChan; jf++) {
	if ((avg[indx+jf]>-9900.0) && (meanRMS>0.0) &&
	    (RMS[indx+jf] > 
	     sqrt(maxRMS[0]*meanRMS*maxRMS[0]*meanRMS + 
		  (avg[indx+jf]*maxRMS[1]*avg[indx+jf]*maxRMS[1]))))
	  avg[indx+jf] = -9999.0;
      } /* end flagging loop */
    } /* end RMS flagging */
    
    /* Get residual RMS */
    sum1 = sum2 = 0.0; cnt = 0.0;
    for (jf=0; jf<numChan; jf++) {
      if (avg[indx+jf]>-9900.0) {
	if (avg[indx+jf]>0.0) {
	  cnt++; 
	  sum1 += avg[indx+jf];
	  sum2 += avg[indx+jf]*avg[indx+jf];
	}
      }
    }
    /* Residual RMS */
    if (cnt>0) residualRMS = sqrt (sum2/cnt + (sum1/cnt)*(sum1/cnt));
    else residualRMS = 1.0e20;
    
    /* Residual flagging */
    for (jf=0; jf<numChan; jf++) {
      if (fabs(avg[indx+jf])>maxRes*sigma[indx+jf]) avg[indx+jf] = -9999.0;
    } /* end loop over Channel */
    jndx += numChan;
    indx += numChan;
  } /* end loop over spectrum */

  /* Flag data with fewer than minGood of time samples */
  /* Find maximum number of time samples */
  maxCount = 0;
  for (jsp=loSp; jsp<=hiSp; jsp++) {
    indx = jsp*numChan;
    for (jf=0; jf<numChan; jf++) {
      maxCount = MAX (count[indx++], maxCount);
    }
  } /* nd loop over spectra */

  flagCount = (ollong)(maxCount*minGood);  /* Level at which to flag */
  if (flagCount>0) {
    for (jsp=loSp; jsp<=hiSp; jsp++) {
     indx = jsp*numChan;
     for (jf=0; jf<numChan; jf++) {
       if ((count[indx]<=flagCount) && (count[indx]>0)) avg[indx] = -9999.0;
       indx++;
     }
    } /* end loop over spectra */
  } /* end if enough time samples to flag */

  /* Indicate completion */
 done:
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadEditFDProcess */

/**
 * Determine robust RMS value of a ofloat array about mean
 * Use center 1-alpha of points, excluding at least one point from each end
 * \param n       Number of points, needs at least 4
 * \param value   Array of values assumed sorted
 * \param mean    Mean value of value
 * \param alpha   0 -> 1 = pure boxcar -> pure MWF (ALPHA of the 
 *                data samples are discarded and the rest averaged). 
 * \return RMS value, fblank if cannot determine
 */
static ofloat FDMedianSigma (olong n, ofloat *value, ofloat mean, ofloat alpha)
{
  ofloat fblank = ObitMagicF();
  ofloat out;
  ofloat sum, beta;
  olong i, i1, i2, count;

  out = fblank;
  if (n<=4) return out;
  if (mean==fblank) return out;

  beta = MAX (0.05, MIN (0.95, 1.0-alpha)) / 2.0; /*  Average around median factor */

  /* Get RMS around the center 1-alpha */
  i1 = MAX (1,   (n/2)-(olong)(beta*n+0.5));
  i2 = MIN (n-1, (n/2)+(olong)(beta*n+0.5));

  if (i2>i1) {
    sum = 0.0;
    count = 0;
    for (i=i1; i<i2; i++) {
      if (value[i]!=fblank) {
	sum += (value[i]-mean)*(value[i]-mean);
	count++;
      }
    }
    if (count>1) out = sqrt(sum / (count-1));
  }
   
  return out;
} /* end FDMedianSigma */

