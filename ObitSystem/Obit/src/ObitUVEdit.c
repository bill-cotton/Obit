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
#include "ObitUVUtil.h"
#include "ObitTableFG.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVEdit.c
 * ObitUVEdit module function definitions.
 * Routines for automatically editing data.
 */

/*----------------- Macroes ---------------------------*/
/*---------------Private structures----------------*/
/* Median editing threaded function argument */
typedef struct {
  /* ObitThread with restart queue */
  ObitThread *thread;
  /* thread number, <0 -> no threading   */
  olong        ithread;
  /* First (1-rel) baseline to process this thread */
  olong        first;
  /* Highest (1-rel) baseline to process this thread  */
  olong        last;
  /* Circular amplitude storage (baseline,correlator,time) */
  ofloat *amps;
  /*  time slice of interest */
  olong itime;  
  /* times of data samples, <-1000 => ignore */
  ofloat* times;   
  /*  Number of baselines */
  olong numBL;
  /* Number of correlations */
  olong numCorr; 
  /* Number of allocated time slices in amps */
  olong numTime; 
  /*  Number of actual time slices */
  olong ntime;  
  /* Controls averaging 
     0 -> 1 = pure boxcar -> pure MWF (alpha of the 
     center data samples are ignored and the rest averaged).
  */
  ofloat alpha;   
  /* [out] Deviations in sigma per baseline/correlator, fblank => not determined */
  ofloat* devs;   
  /* Work array at least the size of ntime */
  ofloat* work;
  /* [out] number of baselines/correlations with valid data */
  ollong number; 
} UVMednEditFuncArg;

/*---------------Private function prototypes----------------*/
/** 
 * Private: Determine minimum clipping levels based on histogram of baseline  
 * RMSes. for TD
 */
static void editHist (olong ncorr, olong numCell, ofloat hisinc, olong *hissig, 
		      ofloat* hiclip);

/** Digest correlator information for TD */
static void digestCorr (ObitUVDesc *inDesc, ofloat *maxRMS, ofloat *maxRMS2, 
			olong *crossBL1, olong *crossBL2, 
			olong *corChan, olong *corIF, olong *corStok);

/** Digest correlator information for TDRMSAvg */
static void digestCorrTDRMSAvg (ObitUVDesc *inDesc, ofloat maxRMSAvg, ofloat *maxRMS2, 
				olong *crossBL1, olong *crossBL2, 
				olong *corChan, olong *corIF, olong *corStok);

/**Private:  Digest  some correlator information for FD */
static void digestCorrFD (ObitUVDesc *inDesc, olong *corChan, olong *corIF, 
			 olong *corStok, olong *corV);

/** Private: FD spectral baseline mask  */
static void EditFDChanMask(olong numChan, olong numIF, olong numPol, olong *chanSel, 
			   gboolean *chanMask);

/** Private: Average for FD, get RMSes  */
static void EditFDAvg(olong numChan, olong numIF, olong numPol, olong numBL, ollong *count, 
		      ofloat *sumA, ofloat *sumA2, ofloat maxAmp, 
		      ofloat maxV, olong *corV);

/** Private: Fit baselines  for FD to averages, subtract */
static void EditFDBaseline(olong numChan, olong numIF, olong numPol, olong numBL, 
			   ollong *count, ofloat *avg, ofloat *RMS, ofloat *sigma,
			   olong widMW, gboolean *chanMask, ObitErr* err);
	    
/** Private: Do editing for FD */
static void EditFDEdit(olong numChan, olong numIF, olong numPol, olong numBL, 
		       ollong *count, ofloat *avg, ofloat *RMS, ofloat *sigma,
		       ofloat *maxRMS, ofloat maxRes, ofloat maxResBL, 
		       gboolean *chanMask, gboolean doMW);

/** Private: Get median value of an array */
static ofloat medianVal (ofloat *array, olong incs, olong n);

/** Private: qsort ofloat comparison */
static int compare_ofloat  (const void* arg1,  const void* arg2);

/** Private: Fit linear slope with mask */
static void EditFDBLFit (ofloat* x, ofloat* y, gboolean *m, olong ndata, 
			 ofloat *a, ofloat *b);

/** Private: Determine deviations from Median */
static ollong MedianDev (ofloat *amps, olong itime, ofloat *times,
			 olong numBL, olong ncorr, olong numTime, olong ntime,
			 ofloat alpha, ofloat *devs, 
			 ofloat *work, olong nThread, UVMednEditFuncArg** args,
			 ObitErr *err);

/** Private: Return single uv vis with buffered read */
static ObitIOCode ReadOne (ObitUV *inUV, gboolean doCalSelect, ofloat** Buffer, 
			   olong *visNo, ObitErr *err);
/** Private: Determine data integration */
static ofloat MedianUVInt (ObitUV *inUV, ObitErr *err);

/** Private: Determine average level for Median */
static ofloat MedianLevel (olong n, ofloat *value, ofloat alpha);

/** Private: Determine sigma for Median */
static ofloat MedianSigma (olong n, ofloat *value, ofloat mean, ofloat alpha);

/** Private: Median flagging */
static ollong MedianFlag (ofloat *devs, ofloat flagSig, 
			  olong numBL, olong numCorr, gboolean allFlag,
			  ofloat time, ofloat timeInt,
			  olong *BLAnt1, olong *BLAnt2, 
			  olong *Chan, olong *IFs, olong *Stoke,
			  ObitTableFG *FGtab, ObitTableFGRow *FGRow, 
			  ObitErr *err);

/** Private: Get array of channel ranges for median flagging */
static olong* medianChan (ObitUVDesc *inDesc, ObitUVDesc *outDesc, olong begChan);

/** Private: Get array of IF ranges for median flagging */
static olong* medianIFs (ObitUVDesc *inDesc, ObitUVDesc *outDesc, olong begIF);

/** Private: Get array of Poln flags for median flagging */
static olong* medianPoln (ObitUVDesc *inDesc, ObitUVDesc *outDesc);
/** Private: Time to String */
static void T2String (ofloat time, gchar *msgBuf);
/** Private: Threaded MedianDev */
static gpointer ThreadMedianDev (gpointer arg);


/*----------------------Public functions---------------------------*/
/**
 * Time-domain editing of UV data - produces FG table
 * Fill flagging table with clipping by RMS values of the real and imaginary
 * parts.  All correlations are clipped on each baseline if the RMS is
 * larger than  the maximum.  The clipping is done independently in
 * each time interval defined by timeAvg. 
 *    The clipping level is given by MIN (A, MAX (B,C)) where:
 * A = sqrt (maxRMS[0]**2 + (avg_amp * maxRMS[1])**2)
 *     and avg_amp is the average amplitude on each baseline.
 * B = median RMS + 3 * sigma of the RMS distribution.
 * C = level corresponding to 3% of the data.
 * Usage of C is controlled by doHiEdit
 *    All data on a given baseline/correlator are flagged if the RMS
 * exceeds the limit.  If a fraction of bad baselines on any correlator
 * exceeds maxBad, then all data to that correlator is flagged.  In
 * addition, if the offending correlator is a parallel hand correlator
 * then any corresponding cross hand correlations are also flagged.
 * Flagging entries are written into FG table flagTab.
 * Control parameters on info member of inUV:
 * \li "flagTab" OBIT_int    (1,1,1) FG table version number [ def. 1]
 * \li "timeAvg" OBIT_float  (1,1,1) Time interval over which to determine 
 *               data to be flagged (min) [def = 1 min.]
 *               NB: this should be at least 2 integrations.
 * \li "maxRMS"  OBIT_float (2,1,1) Maximum RMS allowed, constant plus 
 *               amplitude coefficient. 
 * \li "maxBad"  OBIT_float (1,1,1) Fraction of allowed flagged baselines 
 *               to a poln/channel/IF above which all baselines are flagged.
 *               [default 0.25]
 * \li "doHiEdit" OBIT_boolean (1,1,1) If present and TRUE, flag based 
 *                on statistics [def FALSE]
 * Routine adapted from the AIPSish TDEDIT.FOR/TDEDIT
 * \param inUV     Input uv data to edit. 
 * \param outUV    UV data onto which the FG table is to be attached.
 *                 May be the same as inUV.
 * \param err      Error stack, returns if  error.
 */
void ObitUVEditTD (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  ObitTableFG *outFlag=NULL;
  ObitTableFGRow *row=NULL;
  gboolean doCalSelect, doHiEdit;
  olong i, j, k, jj, kk, firstVis, startVis, suba, iFGRow, ver;
  ollong countAll, countBad;
  olong lastSourceID, curSourceID, lastSubA, lastFQID=-1;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  ObitUVDesc *inDesc;
  ofloat timeAvg, lastTime=-1.0, maxRMS[2], *maxRMS2=NULL, maxBad, cbase, hisicc;
  olong *hissig=NULL, *blLookup=NULL, *crossBL1=NULL, *crossBL2=NULL;
  olong *corChan=NULL, *corIF=NULL, *corStok=NULL;
  olong *BLAnt1=NULL, *BLAnt2=NULL, BIF, BChan;
  olong flagTab, indx, jndx, kndx, nVisPIO, itemp, ant1, ant2;
  olong numCell, ncorr, numAnt, numBL, blindx, hicell;
  gboolean gotOne, done, isbad, *badCor=NULL;
  ofloat *acc=NULL, *corCnt=NULL, *corBad=NULL, *hiClip=NULL, *Buffer;
  ofloat startTime, endTime, curTime, sigma, hisinc, rms2, rms3, ampl2;
  ofloat mxamp2, sec;
  gchar *tname, reason[25];
  struct tm *lp;
  time_t clock;
  gchar *routine = "ObitUVEditTD";

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
  flagTab = 1;
  ObitInfoListGetTest(inUV->info, "flagTab", &type, dim, &flagTab);
  /* Time interval */
  timeAvg = 1.0;  /* default 1 min */
  ObitInfoListGetTest(inUV->info, "timeAvg", &type, dim, &timeAvg);
  if (timeAvg<=(1.0/60.0)) timeAvg = 1.0;
  timeAvg /= 1440.0;  /* convert to days */
  /* RMS clipping parameters, no default */
  ObitInfoListGet(inUV->info, "maxRMS", &type, dim, maxRMS, err);
  mxamp2 = maxRMS[1] * maxRMS[1];
  /* max. fraction bad baselines */
  maxBad = 0.25;           /* default 0.25 */
  ObitInfoListGetTest(inUV->info, "maxBad", &type, dim,  &maxBad);  
  /* Edit by histogram? */
  doHiEdit = FALSE;
  ObitInfoListGetTest(inUV->info, "doHiEdit", &type, dim,  &doHiEdit);  

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

  /* Selection of input? */
  doCalSelect = TRUE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, inUV->name);
  inDesc  = inUV->myDesc;  /* Get descriptor */

  /* Create output FG table */
  tname = g_strconcat ("Flag Table for: ", outUV->name, NULL);
  ver = flagTab;
  outFlag = newObitTableFGValue(tname, (ObitData*)outUV, &ver, OBIT_IO_ReadWrite, 
				err);
  g_free (tname);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
 
  /* Allocate arrays */
  numCell = 800;  /* number of cells in histogram */
  suba    = 1;
  numAnt  = inUV->myDesc->numAnt[suba-1];/* actually highest antenna number */
  numBL   = (numAnt*(numAnt-1))/2;
  ncorr   = inUV->myDesc->ncorr;
  /* acc index = type + corr * (6) + BL * (6*ncorr)  where BL = 0-rel baseline index */
  acc    = g_malloc0 (ncorr * 6 * numBL * sizeof(ofloat));
  /* hissig index  = cell + numCell*corr */
  hissig = g_malloc0 (ncorr * numCell * sizeof(olong));
  hisicc = 0.005;  /* Histogram resolution ??? This seems critical */
  hisinc = hisicc*maxRMS[0]; /* Histogram increment */
  hiClip = g_malloc0 (ncorr * sizeof(ofloat));
  corCnt = g_malloc0 (ncorr * sizeof(ofloat));
  corBad = g_malloc0 (ncorr * sizeof(ofloat));
  badCor = g_malloc0 (ncorr * sizeof(gboolean));
  maxRMS2  = g_malloc0 (ncorr * sizeof(ofloat)); /* Maximum variance per correlator */
  crossBL1 = g_malloc0 (ncorr * sizeof(olong));   /* associated cross baseline 1 */
  crossBL2 = g_malloc0 (ncorr * sizeof(olong));   /* associated cross baseline 2 */
  corChan  = g_malloc0 (ncorr * sizeof(olong));   /* Correlator Channel */
  corIF    = g_malloc0 (ncorr * sizeof(olong));   /* Correlator IF */
  corStok  = g_malloc0 (ncorr * sizeof(olong));   /* Correlator Stokes */

  /* Baseline tables */
  blLookup = g_malloc0(numAnt*sizeof(olong));
  BLAnt1   = g_malloc0 (numBL * sizeof(olong));    /* First antenna of baseline */
  BLAnt2   = g_malloc0 (numBL * sizeof(olong));    /* Second antenna of baseline */
  blLookup[0] = 0;
  k = 0;
  for (j=2; j<=numAnt; j++) {BLAnt1[k]=1; BLAnt2[k]=j; k++;}
  for (i=1; i<numAnt; i++) {
    blLookup[i] = blLookup[i-1] + numAnt-i;
    for (j=i+2; j<=numAnt; j++) {BLAnt1[k]=i+1; BLAnt2[k]=j; k++;}
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
  countBad     = 0;
  for (i=0; i<6*ncorr*numBL; i++) acc[i] = 0.0;

  Buffer = inUV->buffer;

  /* Digest visibility info */
  digestCorr (inDesc, maxRMS, maxRMS2, crossBL1, crossBL2, 
	    corChan, corIF, corStok);

  /* Open output table */
  oretCode = ObitTableFGOpen (outFlag, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  
  /* Create Row */
  row = newObitTableFGRow (outFlag);
  
  /* Attach row to output buffer */
  ObitTableFGSetRow (outFlag, row, err);
  if (err->error) goto cleanup;
  
  /* If there are entries in the table, mark it unsorted */
  if (outFlag->myDesc->nrow>0) 
    {outFlag->myDesc->sort[0]=0; outFlag->myDesc->sort[1]=0;}
  
  /* Initialize solution row */
  row->SourID  = 0; 
  row->SubA    = 0; 
  row->freqID  = 0; 
  row->ants[0] = 0; 
  row->ants[1] = 0; 
  row->TimeRange[0] = -1.0e20; 
  row->TimeRange[1] =  1.0e20; 
  row->ifs[0]    = BIF; 
  row->ifs[1]    = 0; 
  row->chans[0]  = BChan; 
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
  g_snprintf (reason, 25, "EditTD %d/%d/%d %d:%d:%f", 
	      lp->tm_year, lp->tm_mon, lp->tm_mday, 
	      lp->tm_hour, lp->tm_min, sec);
  row->reason    = reason; /* Unique string */
  
  /* Loop over intervals */
  done   = FALSE;
  gotOne = FALSE;
  while (!done) {
    
    /* we're in business - loop through data - one vis per read */
    while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
      if ((!gotOne) || (inUV->myDesc->numVisBuff<=0)) { /* need to read new record? */
	if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
	else iretCode = ObitUVRead (inUV, inUV->buffer, err);
	if (err->error) goto cleanup;
	/*if (iretCode!=OBIT_IO_OK) break;*/
	firstVis = inDesc->firstVis;
      }
      
      /* Are we there yet??? */
      done = (inDesc->firstVis >= inDesc->nvis) || 
	(iretCode==OBIT_IO_EOF);
      if (done && (startTime>0.0)) goto process;


      /* Make sure valid data found */
      if (inUV->myDesc->numVisBuff<=0) continue;
      
      gotOne = FALSE;
      /* Time */
      curTime = Buffer[inDesc->iloct];
      if (inDesc->ilocsu>=0) curSourceID = Buffer[inDesc->ilocsu];
      if (startTime < -1000.0) {  /* Set time window etc. if needed */
	startTime = curTime;
	lastTime = curTime;
	endTime   = startTime +  timeAvg;
	startVis  = firstVis;
	lastSourceID = curSourceID;
      }

      /* Still in current interval/source? */
      if ((curTime<endTime) && (curSourceID == lastSourceID) && 
	  (inDesc->firstVis<=inDesc->nvis) && (iretCode==OBIT_IO_OK)) {
	/* accumulate statistics */
	cbase = Buffer[inUV->myDesc->ilocb]; /* Baseline */
	ant1 = (cbase / 256.0) + 0.001;
	ant2 = (cbase - ant1 * 256) + 0.001;
	lastSubA = (olong)(100.0 * (cbase -  ant1 * 256 - ant2) + 1.5);
	/* Baseline index this assumes a1<a2 always - ignore auto correlations */
	if (ant1!=ant2) {
	  blindx =  blLookup[ant1-1] + ant2-ant1-1;
	  if (inDesc->ilocfq>=0) lastFQID = (olong)(Buffer[inDesc->ilocfq]+0.5);
	  else lastFQID = 0;
	  lastTime = curTime;
	  
	  /* Accumulate
	     (1,*) =  count 
	     (2,*) =  sum r then max rms**2 
	     (3,*) =  sum**2 real then rms**2 
	     (4,*) =  sum imaginary 
	     (5,*) =  sum**2 imaginary then rms**2 
	     (6,*) =  sum amplitude */
	  indx = inDesc->nrparm; /* offset of start of vis data */
	  for (i=0; i<ncorr; i++) { /* loop 120 */
	    if (Buffer[indx+2] > 0.0) {
	      jndx = i*6 + blindx*6*ncorr;
	      acc[jndx]   += 1.0;
	      acc[jndx+1] += Buffer[indx];
	      acc[jndx+2] += Buffer[indx] * Buffer[indx];
	      acc[jndx+3] += Buffer[indx+1];
	      acc[jndx+4] += Buffer[indx+1] * Buffer[indx+1];
	      acc[jndx+5] += sqrt (Buffer[indx]*Buffer[indx] + Buffer[indx+1]*Buffer[indx+1]);
	    } 
	    indx += 3;
	  } /* end loop  L120 over correlations */;
	} /* end only cross correlations */
      } else {  /* process interval */
	
      process:
	/* Now have the next record in the IO Buffer */
	if (iretCode==OBIT_IO_OK) gotOne = TRUE;
	    
	/* Get RMS/collect histogram */
	for (i=0; i<ncorr * numCell; i++)  hissig[i] = 0;
	for (i=0; i<numBL; i++) { /* loop 170 */
	  for (j=0; j<ncorr; j++) { /* loop 160 */
 	    jndx = j*6 + i*6*ncorr;
	    if (acc[jndx] > 0.0) countAll++;  /* count all possibilities */
	    if (acc[jndx] > 1.1) {
	      /* Real part */
	      rms2 = (acc[jndx+2] - ((acc[jndx+1]*acc[jndx+1]) / acc[jndx])) / (acc[jndx]-1.0);
	      rms2 = fabs (rms2);
	      acc[jndx+2] = rms2;
	      /* Imaginary part */
	      rms3 = (acc[jndx+4] - ((acc[jndx+3]*acc[jndx+3]) / acc[jndx])) / (acc[jndx]-1.0);
	      rms3 = fabs (rms3);
	      acc[jndx+4] = rms3;
	      /* Histogram */
	      sigma = sqrt (MAX (rms2, rms3));
	      hicell = sigma / hisinc;
	      hicell = MIN (hicell, numCell-1);
	      hissig[hicell+j*numCell]++;
	    } 
	  } /* end loop  L160: */
	} /* end loop  L170: */

	if (doHiEdit) {
	/* Set histogram clipping levels. */
	  editHist (ncorr, numCell, hisinc, hissig, hiClip);
	} else {
	  /* No histogram flagging  */
	  for (j=0; j<ncorr; j++) { /* loop 180 */
            hiClip[j] = maxRMS[0] * maxRMS[0];
	  } /* end loop  L180: */;
	}

	/* initialize counters */
	for (i=0; i<ncorr; i++) corCnt[i] = corBad[i] = 0;

	/* Find bad baselines. */
	for (i=0; i<numBL; i++) { /* loop 200 */
	  for (j=0; j<ncorr; j++) { /* loop 190 */
 	    jndx = j*6 + i*6*ncorr;
            if (acc[jndx] > 1.1) {
	      /* Real part */
	      rms2 = acc[jndx+2];
	      /* Convert sum to maximum rms**2 */
	      ampl2 = (acc[jndx+5] / acc[jndx]) * (acc[jndx+5] / acc[jndx]);
	      acc[jndx+1] = MIN ((maxRMS2[j] + (mxamp2*ampl2)), hiClip[j]);
	      /* Is this one bad? */
	      isbad = rms2  >  acc[jndx+1];
	      /* Imaginary part */
	      rms2 = acc[jndx+4];
	      /* Convert sum to maximum rms**2 */
	      acc[jndx+3] = MIN ((maxRMS2[j] + (mxamp2*ampl2)), hiClip[j]);
	      /* Is this one bad? */
	      isbad = isbad  ||  (rms2  >  acc[jndx+3]);

	      /* Correlator info */
	      corCnt[j]++;
	      if (isbad) {
		/* Make sure it is flagged. */
		acc[jndx+2] = 1.0e20;
		acc[jndx+1] = 0.0;
		corBad[j]++;
		/* If parallel and bad, kill its crosspolarized relatives */
		if ((crossBL1[j]>0) && (crossBL1[i]<ncorr)) {
		  kndx = crossBL1[j]*6 + i*6*ncorr;
		  acc[kndx+2] = 1.0e20;
		  acc[kndx+1] = 0.0;
		}
		if ((crossBL2[j]>0) && (crossBL2[i]<ncorr)) {
		  kndx = crossBL2[j]*6 + i*6*ncorr;
		  acc[kndx+2] = 1.0e20;
		  acc[kndx+1] = 0.0;
		}
	      } 
	    } else if (acc[jndx] > 0.0) {
	      /* Flag correlators without enough data. */
	      acc[jndx+2] = 1.0e20;
	      acc[jndx+1] = 0.0;
	    }
	  } /* end loop  L190: */
	} /* end loop  L200: */
	
	/* Check for bad correlators */
	for (i=0; i<ncorr; i++) { /* loop 210 */
	  badCor[i] = FALSE;
	  if (corCnt[i] > 1.1) {
	    if ((corBad[i]/corCnt[i])  >  maxBad) {
	      /* Kill correlator... */
	      badCor[i] = TRUE;
	      /* and all its relatives */
	      if ((crossBL1[i]>0) && (crossBL1[i]<ncorr)) 
		badCor[crossBL1[i]] = TRUE;
	      if ((crossBL2[i]>0) && (crossBL2[i]<ncorr)) 
		badCor[crossBL2[i]] = TRUE;
	    }
	  } 
	} /* end loop  L210: */
	
	
	/* Init Flagging table entry */
	row->SourID  = lastSourceID; 
	row->SubA    = lastSubA; 
	row->freqID  = lastFQID; 
	row->TimeRange[0] = startTime;
	row->TimeRange[1] = lastTime;
	
	/* Loop over correlators flagging bad */
	for (i=0; i<ncorr; i++) { /* loop 210 */
	  if (badCor[i]) {
	    row->ants[0]  = 0; 
	    row->ants[1]  = 0; 
	    row->ifs[0]   = BIF + corIF[i] - 1; 
	    row->ifs[1]   = BIF + corIF[i] - 1; 
	    row->chans[0] = BChan + corChan[i]  - 1; 
	    row->chans[1] = BChan + corChan[i] - 1; 
	    row->pFlags[0]=row->pFlags[1]=row->pFlags[2]=row->pFlags[3]=0; 
	    k = abs (corStok[i]);
	    /* bit flag implementation kinda screwy */
	    row->pFlags[0] |= 1<<(k-1);
	    
	    /* Write */
	    iFGRow = -1;
	    oretCode = ObitTableFGWriteRow (outFlag, iFGRow, row, err);
	    if (err->error) goto cleanup;
	  } /* end bad correlator section */
	} /* end loop flagging correlators */
	
	/* Loop over baselines/correlator flagging bad */
	for (i=0; i<numBL; i++) {
	  for (j=0; j<ncorr; j++) {
	    jndx = j*6 + i*6*ncorr;
	    /* Count flagged interval/correlator */
	    if ((acc[jndx]>0.0) && (badCor[j] || (acc[jndx+2] > acc[jndx+1])))
	      countBad++;
	    if ((!badCor[j])) {  /* Don't flag if correlator already flagged */
	      if ((acc[jndx]>0.0) && (acc[jndx+2] > acc[jndx+1])) {
		/* Check for higher number spectral channels in a contigious 
		   block of bad channels and include in same flag */
		jj = j;
		for (kk = j+1; kk<ncorr; kk++) {
		  /* Only interested in same IF, poln */
		  if ((corIF[j]!=corIF[kk]) || (corStok[j]!=corStok[kk])) continue;
		  kndx = kk*6 + i*6*ncorr;
		  /* Higher channel number and to be flagged? */
		  if ((corChan[kk]>corChan[jj]) && 
		      ((acc[kndx]>0.0) && (acc[kndx+2]>acc[kndx+1]))) {
		    jj = kk;         /* This correlator to be included in flag */
		    countBad++;      /* Increment bad count */
		    acc[kndx] = 0.0; /* don't consider again */
		  } else if (corChan[kk]>corChan[jj]) { /* Higher channel number and good? */
		    break;  /* stop looking at higher channel numbers */
		  }
		} /* end loop searching for higher bad channels */
		row->ants[0]   = BLAnt1[i]; 
		row->ants[1]   = BLAnt2[i]; 
		row->ifs[0]    = BIF + corIF[j] - 1; 
		row->ifs[1]    = BIF + corIF[j] - 1; 
		row->chans[0]  = BChan + corChan[j]  - 1; 
		row->chans[1]  = BChan + corChan[jj] - 1; 
		row->pFlags[0]=row->pFlags[1]=row->pFlags[2]=row->pFlags[3]=0; 
		k = abs (corStok[j]);
		/* bit flag implementation kinda screwy */
		row->pFlags[0] |= 1<<(k-1);
		
		/* Write */
		iFGRow = -1;
		oretCode = ObitTableFGWriteRow (outFlag, iFGRow, row, err);
		if (err->error) goto cleanup;
	      }  /* end flag correlator */
	    } /* end correlator not flagged */
	  } /* end loop over correlators */
	} /* end loop over baselines */
	
	/* Are we there yet??? */
	done = (inDesc->firstVis >= inDesc->nvis) || 
	  (iretCode==OBIT_IO_EOF);

	/* Reinitialize things */
	startTime = -1.0e20;
	endTime   =  1.0e20;
	for (i=0; i<6*ncorr*numBL; i++) acc[i] = 0.0;

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

  /* Cleanup */
 cleanup:  
  /* close uv file */
  iretCode = ObitUVClose (inUV, err);
  
  /* Close output table */
  oretCode = ObitTableFGClose (outFlag, err);
  
  /* Deallocate arrays */
  row     = ObitTableFGRowUnref(row);
  outFlag = ObitTableFGUnref(outFlag);
  if (acc)      g_free(acc);
  if (hissig)   g_free(hissig);
  if (hiClip)   g_free(hiClip);
  if (corCnt)   g_free(corCnt);
  if (corBad)   g_free(corBad);
  if (badCor)   g_free(badCor);
  if (blLookup) g_free(blLookup);
  if (maxRMS2)  g_free(maxRMS2);
  if (crossBL1) g_free(crossBL1);
  if (crossBL2) g_free(crossBL2);
  if (corChan)  g_free(corChan);
  if (corIF)    g_free(corIF);
  if (corStok)  g_free(corStok);
  if (BLAnt1)   g_free(BLAnt1);
  if (BLAnt2)   g_free(BLAnt2);
  if (err->error)  Obit_traceback_msg (err, routine, inUV->name);

  /* Give report */
  Obit_log_error(err, OBIT_InfoErr, "EditTD: flag %8.1lf of %8.1lf vis/interval= %5.1lf percent",
		 (odouble)countBad, (odouble)countAll, 100.0*(odouble)countBad/((odouble)countAll));

  return;
} /* end  ObitUVEditTD */

/**
 * Time-domain RMS/Avg  editing of UV data - produces FG table
 * Fill flagging table with clipping by (rms_amp/amp)
 * All correlations are clipped on each baseline if the RMS/Avg is
 * larger than  the maximum.  The clipping is done independently in
 * each time interval defined by timeAvg. 
 *    The clipping level is given by maxRMSAvg
 *    All data on a given baseline/correlator are flagged if the RMS
 * exceeds the limit.  If a fraction of bad baselines on any correlator
 * exceeds maxBad, then all data to that correlator is flagged.  In
 * addition, if the offending correlator is a parallel hand correlator
 * then any corresponding cross hand correlations are also flagged.
 * Flagging entries are written into FG table flagTab.
 * Control parameters on info member of inUV:
 * \li "flagTab" OBIT_int    (1,1,1) FG table version number [ def. 1]
 * \li "timeAvg" OBIT_float  (1,1,1) Time interval over which to determine 
 *               data to be flagged (min) [def = 1 min.]
 *               NB: this should be at least 2 integrations.
 * \li "maxRMSAvg" OBIT_float (1,1,1) Maximum RMS/avg allowed
 * \li "maxBad"    OBIT_float (1,1,1) Fraction of allowed flagged baselines 
 *               to a poln/channel/IF above which all baselines are flagged.
 *               [default 0.25]
 *
 * \param inUV     Input uv data to edit. 
 * \param outUV    UV data onto which the FG table is to be attached.
 *                 May be the same as inUV.
 * \param err      Error stack, returns if  error.
 */
void ObitUVEditTDRMSAvg (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  ObitTableFG *outFlag=NULL;
  ObitTableFGRow *row=NULL;
  gboolean doCalSelect;
  olong i, j, k, jj, kk, firstVis, startVis, suba, iFGRow, ver;
  ollong countAll, countBad;
  olong lastSourceID, curSourceID, lastSubA, lastFQID=-1;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  ObitUVDesc *inDesc;
  ofloat timeAvg, lastTime=-1.0, maxRMSAvg, *maxRMS2=NULL, maxBad, cbase;
  olong *blLookup=NULL, *crossBL1=NULL, *crossBL2=NULL;
  olong *corChan=NULL, *corIF=NULL, *corStok=NULL;
  olong *BLAnt1=NULL, *BLAnt2=NULL, BIF, BChan;
  olong flagTab, indx, jndx, kndx, nVisPIO, itemp, ant1, ant2;
  olong ncorr, numAnt, numBL, blindx;
  gboolean gotOne, done, isbad, *badCor=NULL;
  ofloat *acc=NULL, *corCnt=NULL, *corBad=NULL, *Buffer;
  ofloat startTime, endTime, curTime, rms2, ampl2;
  ofloat sec;

  gchar *tname, reason[25];
  struct tm *lp;
  time_t clock;
  gchar *routine = "ObitUVEditTDRMSAvg";

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
  flagTab = 1;
  ObitInfoListGetTest(inUV->info, "flagTab", &type, dim, &flagTab);
  /* Time interval */
  timeAvg = 1.0;  /* default 1 min */
  ObitInfoListGetTest(inUV->info, "timeAvg", &type, dim, &timeAvg);
  if (timeAvg<=(1.0/60.0)) timeAvg = 1.0;
  timeAvg /= 1440.0;  /* convert to days */
  /* RMS clipping parameters, no default */
  ObitInfoListGet(inUV->info, "maxRMSAvg", &type, dim, &maxRMSAvg, err);
  /* max. fraction bad baselines */
  maxBad = 0.25;           /* default 0.25 */
  ObitInfoListGetTest(inUV->info, "maxBad", &type, dim,  &maxBad);  
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

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

  /* Selection of input? */
  doCalSelect = TRUE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, inUV->name);
  inDesc  = inUV->myDesc;  /* Get descriptor */

  /* Create output FG table */
  tname = g_strconcat ("Flag Table for: ", outUV->name, NULL);
  ver = flagTab;
  outFlag = newObitTableFGValue(tname, (ObitData*)outUV, &ver, OBIT_IO_ReadWrite, 
				err);
  g_free (tname);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
 
  /* Allocate arrays */
  suba    = 1;
  numAnt  = inUV->myDesc->numAnt[suba-1];/* actually highest antenna number */
  numBL   = (numAnt*(numAnt-1))/2;
  ncorr   = inUV->myDesc->ncorr;
  /* acc index = type + corr * (6) + BL * (3*ncorr)  where BL = 0-rel baseline index */
  acc    = g_malloc0 (ncorr * 3 * numBL * sizeof(ofloat));
  corCnt = g_malloc0 (ncorr * sizeof(ofloat));
  corBad = g_malloc0 (ncorr * sizeof(ofloat));
  badCor = g_malloc0 (ncorr * sizeof(gboolean));
  maxRMS2  = g_malloc0 (ncorr * sizeof(ofloat)); /* Maximum variance per correlator */
  crossBL1 = g_malloc0 (ncorr * sizeof(olong));   /* associated cross baseline 1 */
  crossBL2 = g_malloc0 (ncorr * sizeof(olong));   /* associated cross baseline 2 */
  corChan  = g_malloc0 (ncorr * sizeof(olong));   /* Correlator Channel */
  corIF    = g_malloc0 (ncorr * sizeof(olong));   /* Correlator IF */
  corStok  = g_malloc0 (ncorr * sizeof(olong));   /* Correlator Stokes */

  /* Baseline tables */
  blLookup = g_malloc0(numAnt*sizeof(olong));
  BLAnt1   = g_malloc0 (numBL * sizeof(olong));    /* First antenna of baseline */
  BLAnt2   = g_malloc0 (numBL * sizeof(olong));    /* Second antenna of baseline */
  blLookup[0] = 0;
  k = 0;
  for (j=2; j<=numAnt; j++) {BLAnt1[k]=1; BLAnt2[k]=j; k++;}
  for (i=1; i<numAnt; i++) {
    blLookup[i] = blLookup[i-1] + numAnt-i;
    for (j=i+2; j<=numAnt; j++) {BLAnt1[k]=i+1; BLAnt2[k]=j; k++;}
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
  countBad     = 0;
  for (i=0; i<3*ncorr*numBL; i++) acc[i] = 0.0;

  Buffer = inUV->buffer;

  /* Digest visibility info */
  digestCorrTDRMSAvg (inDesc, maxRMSAvg, maxRMS2, crossBL1, crossBL2, 
		      corChan, corIF, corStok);

  /* Open output table */
  oretCode = ObitTableFGOpen (outFlag, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  
  /* Create Row */
  row = newObitTableFGRow (outFlag);
  
  /* Attach row to output buffer */
  ObitTableFGSetRow (outFlag, row, err);
  if (err->error) goto cleanup;
  
  /* If there are entries in the table, mark it unsorted */
  if (outFlag->myDesc->nrow>0) 
    {outFlag->myDesc->sort[0]=0; outFlag->myDesc->sort[1]=0;}
  
  /* Initialize solution row */
  row->SourID  = 0; 
  row->SubA    = 0; 
  row->freqID  = 0; 
  row->ants[0] = 0; 
  row->ants[1] = 0; 
  row->TimeRange[0] = -1.0e20; 
  row->TimeRange[1] =  1.0e20; 
  row->ifs[0]    = BIF; 
  row->ifs[1]    = 0; 
  row->chans[0]  = BChan; 
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
  g_snprintf (reason, 25, "EditTDRMS %d/%d/%d %d:%d:%f", 
	      lp->tm_year, lp->tm_mon, lp->tm_mday, 
	      lp->tm_hour, lp->tm_min, sec);
  row->reason    = reason; /* Unique string */
  
  /* Loop over intervals */
  done   = FALSE;
  gotOne = FALSE;
  while (!done) {
    
    /* we're in business - loop through data - one vis per read */
    while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
      if ((!gotOne) || (inUV->myDesc->numVisBuff<=0)) { /* need to read new record? */
	if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
	else iretCode = ObitUVRead (inUV, inUV->buffer, err);
	if (err->error) goto cleanup;
	/*if (iretCode!=OBIT_IO_OK) break;*/
	firstVis = inDesc->firstVis;
      }
      
      /* Are we there yet??? */
      done = (inDesc->firstVis >= inDesc->nvis) || 
	(iretCode==OBIT_IO_EOF);
      if (done && (startTime>0.0)) goto process;


      /* Make sure valid data found */
      if (inUV->myDesc->numVisBuff<=0) continue;
      
      gotOne = FALSE;
      /* Time */
      curTime = Buffer[inDesc->iloct];
      if (inDesc->ilocsu>=0) curSourceID = Buffer[inDesc->ilocsu];
      if (startTime < -1000.0) {  /* Set time window etc. if needed */
	startTime = curTime;
	lastTime = curTime;
	endTime   = startTime +  timeAvg;
	startVis  = firstVis;
	lastSourceID = curSourceID;
      }

      /* Still in current interval/source? */
      if ((curTime<endTime) && (curSourceID == lastSourceID) && 
	  (inDesc->firstVis<=inDesc->nvis) && (iretCode==OBIT_IO_OK)) {
	/* accumulate statistics */
	cbase = Buffer[inUV->myDesc->ilocb]; /* Baseline */
	ant1 = (cbase / 256.0) + 0.001;
	ant2 = (cbase - ant1 * 256) + 0.001;
	lastSubA = (olong)(100.0 * (cbase -  ant1 * 256 - ant2) + 1.5);
	/* Baseline index this assumes a1<a2 always - ignore auto correlations */
	if (ant1!=ant2) {
	  blindx =  blLookup[ant1-1] + ant2-ant1-1;
	  if (inDesc->ilocfq>=0) lastFQID = (olong)(Buffer[inDesc->ilocfq]+0.5);
	  else lastFQID = 0;
	  lastTime = curTime;
	  
	  /* Accumulate
	     (1,*) =  count 
	     (2,*) =  sum amp then max limit
	     (3,*) =  sum amp**2 then rms^2/avg^2 */
	  indx = inDesc->nrparm; /* offset of start of vis data */
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
	    
	/* Get Amp + RMS */
	for (i=0; i<numBL; i++) { /* loop 170 */
	  for (j=0; j<ncorr; j++) { /* loop 160 */
 	    jndx = j*3 + i*3*ncorr;
	    if (acc[jndx] > 0.0) countAll++;  /* count all possibilities */
	    if (acc[jndx] > 1.1) {
	      acc[jndx+1] /= acc[jndx];  /* Average amp */
	      acc[jndx+2] /= acc[jndx];  /* Avg sum amp*amp */
	      /* RMS Amplitude */
	      rms2 = (acc[jndx+2] - acc[jndx+1]*acc[jndx+1]);
	      rms2 = fabs (rms2);
	      acc[jndx+2] = rms2;
	    } 
	  } /* end loop  L160: */
	} /* end loop  L170: */

	/* initialize counters */
	for (i=0; i<ncorr; i++) corCnt[i] = corBad[i] = 0;

	/* Find bad baselines. */
	for (i=0; i<numBL; i++) { /* loop 200 */
	  for (j=0; j<ncorr; j++) { /* loop 190 */
 	    jndx = j*3 + i*3*ncorr;
            if (acc[jndx] > 1.1) {

	      /* Get rms**2 */
	      rms2 = acc[jndx+2];
	      /* Average ampl**2 */
	      ampl2 = acc[jndx+1]*acc[jndx+1];
	      acc[jndx+1] = maxRMS2[j];
	      acc[jndx+2] = rms2 / MAX(ampl2, 1.0e-20);

	      /* Is this one bad? */
	      isbad = acc[jndx+2]  > acc[jndx+1]  ;

	      /* Correlator info */
	      corCnt[j]++;
	      if (isbad) {
		/* Make sure it is flagged. */
		acc[jndx+2] = 1.0e20;
		acc[jndx+1] = 0.0;
		corBad[j]++;
		/* If parallel and bad, kill its cross-polarized relatives */
		if ((crossBL1[j]>0) && (crossBL1[i]<ncorr)) {
		  kndx = crossBL1[j]*3 + i*3*ncorr;
		  acc[kndx+2] = 1.0e20;
		  acc[kndx+1] = 0.0;
		}
		if ((crossBL2[j]>0) && (crossBL2[i]<ncorr)) {
		  kndx = crossBL2[j]*3 + i*3*ncorr;
		  acc[kndx+2] = 1.0e20;
		  acc[kndx+1] = 0.0;
		}
	      } 
	    } else if (acc[jndx] > 0.0) {
	      /* Flag correlators without enough data. */
	      acc[jndx+2] = 1.0e20;
	      acc[jndx+1] = 0.0;
	    }
	  } /* end loop  L190: */
	} /* end loop  L200: */
	
	/* Check for bad correlators */
	for (i=0; i<ncorr; i++) { /* loop 210 */
	  badCor[i] = FALSE;
	  if (corCnt[i] > 1.1) {
	    if ((corBad[i]/corCnt[i])  >  maxBad) {
	      /* Kill correlator... */
	      badCor[i] = TRUE;
	      /* and all its relatives */
	      if ((crossBL1[i]>0) && (crossBL1[i]<ncorr)) 
		badCor[crossBL1[i]] = TRUE;
	      if ((crossBL2[i]>0) && (crossBL2[i]<ncorr)) 
		badCor[crossBL2[i]] = TRUE;
	    }
	  } 
	} /* end loop  L210: */
	
	
	/* Init Flagging table entry */
	row->SourID  = lastSourceID; 
	row->SubA    = lastSubA; 
	row->freqID  = lastFQID; 
	row->TimeRange[0] = startTime;
	row->TimeRange[1] = lastTime;
	
	/* Loop over correlators flagging bad */
	for (i=0; i<ncorr; i++) { /* loop 210 */
	  if (badCor[i]) {
	    row->ants[0]  = 0; 
	    row->ants[1]  = 0; 
	    row->ifs[0]   = BIF + corIF[i] - 1; 
	    row->ifs[1]   = BIF + corIF[i] - 1; 
	    row->chans[0] = BChan + corChan[i]  - 1; 
	    row->chans[1] = BChan + corChan[i] - 1; 
	    row->pFlags[0]= row->pFlags[1] = row->pFlags[2] = row->pFlags[3] = 0; 
	    k = abs (corStok[i]);
	    /* bit flag implementation kinda screwy */
	    row->pFlags[0] |= 1<<(k-1);
	    
	    /* Write */
	    iFGRow = -1;
	    oretCode = ObitTableFGWriteRow (outFlag, iFGRow, row, err);
	    if (err->error) goto cleanup;
	  } /* end bad correlator section */
	} /* end loop flagging correlators */
	
	/* Loop over baselines/correlator flagging bad */
	for (i=0; i<numBL; i++) {
	  for (j=0; j<ncorr; j++) {
	    jndx = j*3 + i*3*ncorr;
	    /* Count flagged interval/correlator */
	    if ((acc[jndx]>1.1) && (badCor[j] || (acc[jndx+2] > acc[jndx+1])))
	      countBad++;
	    if ((!badCor[j])) {  /* Don't flag if correlator already flagged */
	      if ((acc[jndx]>1.1) && (acc[jndx+2] > acc[jndx+1])) {
		/* Check for higher number spectral channels in a contigious 
		   block of bad channels and include in same flag */
		jj = j;
		for (kk = j+1; kk<ncorr; kk++) {
		  /* Only interested in same IF, poln */
		  if ((corIF[j]!=corIF[kk]) || (corStok[j]!=corStok[kk])) continue;
		  kndx = kk*3 + i*3*ncorr;
		  /* Higher channel number and to be flagged? */
		  if ((corChan[kk]>corChan[jj]) && 
		      ((acc[kndx]>1.1) && (acc[kndx+2]>acc[kndx+1]))) {
		    jj = kk;         /* This correlator to be included in flag */
		    countBad++;      /* Increment bad count */
		    acc[kndx] = 0.0; /* don't consider again */
		  } else if (corChan[kk]>corChan[jj]) { /* Higher channel number and good? */
		    break;  /* stop looking at higher channel numbers */
		  }
		} /* end loop searching for higher bad channels */
		row->ants[0]   = BLAnt1[i]; 
		row->ants[1]   = BLAnt2[i]; 
		row->ifs[0]    = BIF + corIF[j] - 1; 
		row->ifs[1]    = BIF + corIF[j] - 1; 
		row->chans[0]  = BChan + corChan[j]  - 1; 
		row->chans[1]  = BChan + corChan[jj] - 1; 
		row->pFlags[0]=row->pFlags[1]=row->pFlags[2]=row->pFlags[3]=0; 
		k = abs (corStok[j]);
		/* bit flag implementation kinda screwy */
		row->pFlags[0] |= 1<<(k-1);
		
		/* Write */
		iFGRow = -1;
		oretCode = ObitTableFGWriteRow (outFlag, iFGRow, row, err);
		if (err->error) goto cleanup;
	      }  /* end flag correlator */
	    } /* end correlator not flagged */
	  } /* end loop over correlators */
	} /* end loop over baselines */
	
	/* Are we there yet??? */
	done = (inDesc->firstVis >= inDesc->nvis) || 
	  (iretCode==OBIT_IO_EOF);

	/* Reinitialize things */
	startTime = -1.0e20;
	endTime   =  1.0e20;
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

  /* Cleanup */
 cleanup:  
  /* close uv file */
  iretCode = ObitUVClose (inUV, err);
  
  /* Close output table */
  oretCode = ObitTableFGClose (outFlag, err);
  
  /* Deallocate arrays */
  row     = ObitTableFGRowUnref(row);
  outFlag = ObitTableFGUnref(outFlag);
  if (acc)      g_free(acc);
  if (corCnt)   g_free(corCnt);
  if (corBad)   g_free(corBad);
  if (badCor)   g_free(badCor);
  if (blLookup) g_free(blLookup);
  if (maxRMS2)  g_free(maxRMS2);
  if (crossBL1) g_free(crossBL1);
  if (crossBL2) g_free(crossBL2);
  if (corChan)  g_free(corChan);
  if (corIF)    g_free(corIF);
  if (corStok)  g_free(corStok);
  if (BLAnt1)   g_free(BLAnt1);
  if (BLAnt2)   g_free(BLAnt2);
  if (err->error)  Obit_traceback_msg (err, routine, inUV->name);

  /* Give report */
  Obit_log_error(err, OBIT_InfoErr, "%s: flag %8.1lf of %8.1lf vis/interval= %5.1lf percent",
		 routine, (odouble)countBad, (odouble)countAll, 
		 100.0*(odouble)countBad/((odouble)countAll));

  return;
} /* end  ObitUVEditTDRMSAvg */

/**
 * Time-domain vector RMS/Avg editing of UV data - produces FG table
 * Fill flagging table with clipping by (rms_real**2 +rms_imag**2)/ampl**2
 * All correlations are clipped on each baseline if the RMS/Avg is
 * larger than  the maximum.  The clipping is done independently in
 * each time interval defined by timeAvg. 
 *    The clipping level is given by maxRMSAvg
 *    All data on a given baseline/correlator are flagged if the RMS
 * exceeds the limit.  If a fraction of bad baselines on any correlator
 * exceeds maxBad, then all data to that correlator is flagged.  In
 * addition, if the offending correlator is a parallel hand correlator
 * then any corresponding cross hand correlations are also flagged.
 * Flagging entries are written into FG table flagTab.
 * Control parameters on info member of inUV:
 * \li "flagTab" OBIT_int    (1,1,1) FG table version number [ def. 1]
 * \li "timeAvg" OBIT_float  (1,1,1) Time interval over which to determine 
 *               data to be flagged (min) [def = 1 min.]
 *               NB: this should be at least 2 integrations.
 * \li "maxRMSAvg" OBIT_float (1,1,1) Maximum RMS/avg allowed
 * \li "maxBad"    OBIT_float (1,1,1) Fraction of allowed flagged baselines 
 *               to a poln/channel/IF above which all baselines are flagged.
 *               [default 0.25]
 *
 * \param inUV     Input uv data to edit. 
 * \param outUV    UV data onto which the FG table is to be attached.
 *                 May be the same as inUV.
 * \param err      Error stack, returns if  error.
 */
void ObitUVEditTDRMSAvgVec (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  ObitTableFG *outFlag=NULL;
  ObitTableFGRow *row=NULL;
  gboolean doCalSelect;
  olong i, j, k, jj, kk, firstVis, startVis, suba, iFGRow, ver;
  ollong countAll, countBad;
  olong lastSourceID, curSourceID, lastSubA, lastFQID=-1;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  ObitUVDesc *inDesc;
  ofloat timeAvg, lastTime=-1.0, maxRMSAvg, *maxRMS2=NULL, maxBad, cbase;
  olong *blLookup=NULL, *crossBL1=NULL, *crossBL2=NULL;
  olong *corChan=NULL, *corIF=NULL, *corStok=NULL;
  olong *BLAnt1=NULL, *BLAnt2=NULL, BIF, BChan;
  olong flagTab, indx, jndx, kndx, nVisPIO, itemp, ant1, ant2;
  olong ncorr, numAnt, numBL, blindx;
  gboolean gotOne, done, isbad, *badCor=NULL;
  ofloat *acc=NULL, *corCnt=NULL, *corBad=NULL, *Buffer;
  ofloat startTime, endTime, curTime, rms2, rms3, ampl2;
  ofloat sec;

  gchar *tname, reason[25];
  struct tm *lp;
  time_t clock;
  gchar *routine = "ObitUVEditTDRMSAvgVec";

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
  flagTab = 1;
  ObitInfoListGetTest(inUV->info, "flagTab", &type, dim, &flagTab);
  /* Time interval */
  timeAvg = 1.0;  /* default 1 min */
  ObitInfoListGetTest(inUV->info, "timeAvg", &type, dim, &timeAvg);
  if (timeAvg<=(1.0/60.0)) timeAvg = 1.0;
  timeAvg /= 1440.0;  /* convert to days */
  /* RMS clipping parameters, no default */
  ObitInfoListGet(inUV->info, "maxRMSAvg", &type, dim, &maxRMSAvg, err);
  /* max. fraction bad baselines */
  maxBad = 0.25;           /* default 0.25 */
  ObitInfoListGetTest(inUV->info, "maxBad", &type, dim,  &maxBad);  
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

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

  /* Selection of input? */
  doCalSelect = TRUE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, inUV->name);
  inDesc  = inUV->myDesc;  /* Get descriptor */

  /* Create output FG table */
  tname = g_strconcat ("Flag Table for: ", outUV->name, NULL);
  ver = flagTab;
  outFlag = newObitTableFGValue(tname, (ObitData*)outUV, &ver, OBIT_IO_ReadWrite, 
				err);
  g_free (tname);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
 
  /* Allocate arrays */
  suba    = 1;
  numAnt  = inUV->myDesc->numAnt[suba-1];/* actually highest antenna number */
  numBL   = (numAnt*(numAnt-1))/2;
  ncorr   = inUV->myDesc->ncorr;
  /* acc index = type + corr * (6) + BL * (5*ncorr)  where BL = 0-rel baseline index */
  acc    = g_malloc0 (ncorr * 5 * numBL * sizeof(ofloat));
  corCnt = g_malloc0 (ncorr * sizeof(ofloat));
  corBad = g_malloc0 (ncorr * sizeof(ofloat));
  badCor = g_malloc0 (ncorr * sizeof(gboolean));
  maxRMS2  = g_malloc0 (ncorr * sizeof(ofloat)); /* Maximum variance per correlator */
  crossBL1 = g_malloc0 (ncorr * sizeof(olong));   /* associated cross baseline 1 */
  crossBL2 = g_malloc0 (ncorr * sizeof(olong));   /* associated cross baseline 2 */
  corChan  = g_malloc0 (ncorr * sizeof(olong));   /* Correlator Channel */
  corIF    = g_malloc0 (ncorr * sizeof(olong));   /* Correlator IF */
  corStok  = g_malloc0 (ncorr * sizeof(olong));   /* Correlator Stokes */

  /* Baseline tables */
  blLookup = g_malloc0(numAnt*sizeof(olong));
  BLAnt1   = g_malloc0 (numBL * sizeof(olong));    /* First antenna of baseline */
  BLAnt2   = g_malloc0 (numBL * sizeof(olong));    /* Second antenna of baseline */
  blLookup[0] = 0;
  k = 0;
  for (j=2; j<=numAnt; j++) {BLAnt1[k]=1; BLAnt2[k]=j; k++;}
  for (i=1; i<numAnt; i++) {
    blLookup[i] = blLookup[i-1] + numAnt-i;
    for (j=i+2; j<=numAnt; j++) {BLAnt1[k]=i+1; BLAnt2[k]=j; k++;}
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
  countBad     = 0;
  for (i=0; i<5*ncorr*numBL; i++) acc[i] = 0.0;

  Buffer = inUV->buffer;

  /* Digest visibility info */
  digestCorrTDRMSAvg (inDesc, maxRMSAvg, maxRMS2, crossBL1, crossBL2, 
		      corChan, corIF, corStok);

  /* Open output table */
  oretCode = ObitTableFGOpen (outFlag, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  
  /* Create Row */
  row = newObitTableFGRow (outFlag);
  
  /* Attach row to output buffer */
  ObitTableFGSetRow (outFlag, row, err);
  if (err->error) goto cleanup;
  
  /* If there are entries in the table, mark it unsorted */
  if (outFlag->myDesc->nrow>0) 
    {outFlag->myDesc->sort[0]=0; outFlag->myDesc->sort[1]=0;}
  
  /* Initialize solution row */
  row->SourID  = 0; 
  row->SubA    = 0; 
  row->freqID  = 0; 
  row->ants[0] = 0; 
  row->ants[1] = 0; 
  row->TimeRange[0] = -1.0e20; 
  row->TimeRange[1] =  1.0e20; 
  row->ifs[0]    = BIF; 
  row->ifs[1]    = 0; 
  row->chans[0]  = BChan; 
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
  g_snprintf (reason, 25, "EditTDRMS %d/%d/%d %d:%d:%f", 
	      lp->tm_year, lp->tm_mon, lp->tm_mday, 
	      lp->tm_hour, lp->tm_min, sec);
  row->reason    = reason; /* Unique string */
  
  /* Loop over intervals */
  done   = FALSE;
  gotOne = FALSE;
  while (!done) {
    
    /* we're in business - loop through data - one vis per read */
    while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
      if ((!gotOne) || (inUV->myDesc->numVisBuff<=0)) { /* need to read new record? */
	if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
	else iretCode = ObitUVRead (inUV, inUV->buffer, err);
	if (err->error) goto cleanup;
	/*if (iretCode!=OBIT_IO_OK) break;*/
	firstVis = inDesc->firstVis;
      }
      
      /* Are we there yet??? */
      done = (inDesc->firstVis >= inDesc->nvis) || 
	(iretCode==OBIT_IO_EOF);
      if (done && (startTime>0.0)) goto process;


      /* Make sure valid data found */
      if (inUV->myDesc->numVisBuff<=0) continue;
      
      gotOne = FALSE;
      /* Time */
      curTime = Buffer[inDesc->iloct];
      if (inDesc->ilocsu>=0) curSourceID = Buffer[inDesc->ilocsu];
      if (startTime < -1000.0) {  /* Set time window etc. if needed */
	startTime = curTime;
	lastTime = curTime;
	endTime   = startTime +  timeAvg;
	startVis  = firstVis;
	lastSourceID = curSourceID;
      }

      /* Still in current interval/source? */
      if ((curTime<endTime) && (curSourceID == lastSourceID) && 
	  (inDesc->firstVis<=inDesc->nvis) && (iretCode==OBIT_IO_OK)) {
	/* accumulate statistics */
	cbase = Buffer[inUV->myDesc->ilocb]; /* Baseline */
	ant1 = (cbase / 256.0) + 0.001;
	ant2 = (cbase - ant1 * 256) + 0.001;
	lastSubA = (olong)(100.0 * (cbase -  ant1 * 256 - ant2) + 1.5);
	/* Baseline index this assumes a1<a2 always - ignore auto correlations */
	if (ant1!=ant2) {
	  blindx =  blLookup[ant1-1] + ant2-ant1-1;
	  if (inDesc->ilocfq>=0) lastFQID = (olong)(Buffer[inDesc->ilocfq]+0.5);
	  else lastFQID = 0;
	  lastTime = curTime;
	  
	  /* Accumulate
	     (1,*) =  count 
	     (2,*) =  sum r then max limit
	     (3,*) =  sum**2 real then rms**2 /avg**2
	     (4,*) =  sum imaginary 
	     (5,*) =  sum**2 imaginary  */
	  indx = inDesc->nrparm; /* offset of start of vis data */
	  for (i=0; i<ncorr; i++) { /* loop 120 */
	    if (Buffer[indx+2] > 0.0) {
	      jndx = i*5 + blindx*5*ncorr;
	      acc[jndx]   += 1.0;
	      acc[jndx+1] += Buffer[indx];
	      acc[jndx+2] += Buffer[indx] * Buffer[indx];
	      acc[jndx+3] += Buffer[indx+1];
	      acc[jndx+4] += Buffer[indx+1] * Buffer[indx+1];
	    } 
	    indx += 3;
	  } /* end loop  L120 over correlations */;
	} /* end only cross correlations */
      } else {  /* process interval */
	
      process:
	/* Now have the next record in the IO Buffer */
	if (iretCode==OBIT_IO_OK) gotOne = TRUE;
	    
	/* Get RMS */
	for (i=0; i<numBL; i++) { /* loop 170 */
	  for (j=0; j<ncorr; j++) { /* loop 160 */
 	    jndx = j*5 + i*5*ncorr;
	    if (acc[jndx] > 0.0) countAll++;  /* count all possibilities */
	    if (acc[jndx] > 1.1) {
	      /* Real part */
	      rms2 = (acc[jndx+2] - ((acc[jndx+1]*acc[jndx+1]) / acc[jndx])) / (acc[jndx]-1.0);
	      rms2 = fabs (rms2);
	      acc[jndx+2] = rms2;
	      /* Imaginary part */
	      rms3 = (acc[jndx+4] - ((acc[jndx+3]*acc[jndx+3]) / acc[jndx])) / (acc[jndx]-1.0);
	      rms3 = fabs (rms3);
	      acc[jndx+4] = rms3;
	    } 
	  } /* end loop  L160: */
	} /* end loop  L170: */

	/* initialize counters */
	for (i=0; i<ncorr; i++) corCnt[i] = corBad[i] = 0;

	/* Find bad baselines. */
	for (i=0; i<numBL; i++) { /* loop 200 */
	  for (j=0; j<ncorr; j++) { /* loop 190 */
 	    jndx = j*5 + i*5*ncorr;
            if (acc[jndx] > 1.1) {
	      /* Get rms */
	      rms2 = acc[jndx+2] + acc[jndx+4];
	      /* Convert sum to maximum rms**2 */
	      ampl2 = (acc[jndx+1]*acc[jndx+1] + acc[jndx+3]*acc[jndx+3]) / (acc[jndx]*acc[jndx]);
	      acc[jndx+1] = maxRMS2[j];
	      acc[jndx+2] = rms2/MAX(ampl2, 1.0e-20);

	      /* Is this one bad? */
	      isbad = acc[jndx+2]  > acc[jndx+1]  ;

	      /* Correlator info */
	      corCnt[j]++;
	      if (isbad) {
		/* Make sure it is flagged. */
		acc[jndx+2] = 1.0e20;
		acc[jndx+1] = 0.0;
		corBad[j]++;
		/* If parallel and bad, kill its crosspolarized relatives */
		if ((crossBL1[j]>0) && (crossBL1[i]<ncorr)) {
		  kndx = crossBL1[j]*5 + i*5*ncorr;
		  acc[kndx+2] = 1.0e20;
		  acc[kndx+1] = 0.0;
		}
		if ((crossBL2[j]>0) && (crossBL2[i]<ncorr)) {
		  kndx = crossBL2[j]*5 + i*5*ncorr;
		  acc[kndx+2] = 1.0e20;
		  acc[kndx+1] = 0.0;
		}
	      } 
	    } else if (acc[jndx] > 0.0) {
	      /* Flag correlators without enough data. */
	      acc[jndx+2] = 1.0e20;
	      acc[jndx+1] = 0.0;
	    }
	  } /* end loop  L190: */
	} /* end loop  L200: */
	
	/* Check for bad correlators */
	for (i=0; i<ncorr; i++) { /* loop 210 */
	  badCor[i] = FALSE;
	  if (corCnt[i] > 1.1) {
	    if ((corBad[i]/corCnt[i])  >  maxBad) {
	      /* Kill correlator... */
	      badCor[i] = TRUE;
	      /* and all its relatives */
	      if ((crossBL1[i]>0) && (crossBL1[i]<ncorr)) 
		badCor[crossBL1[i]] = TRUE;
	      if ((crossBL2[i]>0) && (crossBL2[i]<ncorr)) 
		badCor[crossBL2[i]] = TRUE;
	    }
	  } 
	} /* end loop  L210: */
	
	
	/* Init Flagging table entry */
	row->SourID  = lastSourceID; 
	row->SubA    = lastSubA; 
	row->freqID  = lastFQID; 
	row->TimeRange[0] = startTime;
	row->TimeRange[1] = lastTime;
	
	/* Loop over correlators flagging bad */
	for (i=0; i<ncorr; i++) { /* loop 210 */
	  if (badCor[i]) {
	    row->ants[0]  = 0; 
	    row->ants[1]  = 0; 
	    row->ifs[0]   = BIF + corIF[i] - 1; 
	    row->ifs[1]   = BIF + corIF[i] - 1; 
	    row->chans[0] = BChan + corChan[i]  - 1; 
	    row->chans[1] = BChan + corChan[i] - 1; 
	    row->pFlags[0]=row->pFlags[1]=row->pFlags[2]=row->pFlags[3]=0; 
	    k = abs (corStok[i]);
	    /* bit flag implementation kinda screwy */
	    row->pFlags[0] |= 1<<(k-1);
	    
	    /* Write */
	    iFGRow = -1;
	    oretCode = ObitTableFGWriteRow (outFlag, iFGRow, row, err);
	    if (err->error) goto cleanup;
	  } /* end bad correlator section */
	} /* end loop flagging correlators */
	
	/* Loop over baselines/correlator flagging bad */
	for (i=0; i<numBL; i++) {
	  for (j=0; j<ncorr; j++) {
	    jndx = j*5 + i*5*ncorr;
	    /* Count flagged interval/correlator */
	    if ((acc[jndx]>1.1) && (badCor[j] || (acc[jndx+2] > acc[jndx+1])))
	      countBad++;
	    if ((!badCor[j])) {  /* Don't flag if correlator already flagged */
	      if ((acc[jndx]>1.1) && (acc[jndx+2] > acc[jndx+1])) {
		/* Check for higher number spectral channels in a contigious 
		   block of bad channels and include in same flag */
		jj = j;
		for (kk = j+1; kk<ncorr; kk++) {
		  /* Only interested in same IF, poln */
		  if ((corIF[j]!=corIF[kk]) || (corStok[j]!=corStok[kk])) continue;
		  kndx = kk*5 + i*5*ncorr;
		  /* Higher channel number and to be flagged? */
		  if ((corChan[kk]>corChan[jj]) && 
		      ((acc[kndx]>1.1) && (acc[kndx+2]>acc[kndx+1]))) {
		    jj = kk;         /* This correlator to be included in flag */
		    countBad++;      /* Increment bad count */
		    acc[kndx] = 0.0; /* don't consider again */
		  } else if (corChan[kk]>corChan[jj]) { /* Higher channel number and good? */
		    break;  /* stop looking at higher channel numbers */
		  }
		} /* end loop searching for higher bad channels */
		row->ants[0]   = BLAnt1[i]; 
		row->ants[1]   = BLAnt2[i]; 
		row->ifs[0]    = BIF + corIF[j] - 1; 
		row->ifs[1]    = BIF + corIF[j] - 1; 
		row->chans[0]  = BChan + corChan[j]  - 1; 
		row->chans[1]  = BChan + corChan[jj] - 1; 
		row->pFlags[0]=row->pFlags[1]=row->pFlags[2]=row->pFlags[3]=0; 
		k = abs (corStok[j]);
		/* bit flag implementation kinda screwy */
		row->pFlags[0] |= 1<<(k-1);
		
		/* Write */
		iFGRow = -1;
		oretCode = ObitTableFGWriteRow (outFlag, iFGRow, row, err);
		if (err->error) goto cleanup;
	      }  /* end flag correlator */
	    } /* end correlator not flagged */
	  } /* end loop over correlators */
	} /* end loop over baselines */
	
	/* Are we there yet??? */
	done = (inDesc->firstVis >= inDesc->nvis) || 
	  (iretCode==OBIT_IO_EOF);

	/* Reinitialize things */
	startTime = -1.0e20;
	endTime   =  1.0e20;
	for (i=0; i<5*ncorr*numBL; i++) acc[i] = 0.0;

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

  /* Cleanup */
 cleanup:  
  /* close uv file */
  iretCode = ObitUVClose (inUV, err);
  
  /* Close output table */
  oretCode = ObitTableFGClose (outFlag, err);
  
  /* Deallocate arrays */
  row     = ObitTableFGRowUnref(row);
  outFlag = ObitTableFGUnref(outFlag);
  if (acc)      g_free(acc);
  if (corCnt)   g_free(corCnt);
  if (corBad)   g_free(corBad);
  if (badCor)   g_free(badCor);
  if (blLookup) g_free(blLookup);
  if (maxRMS2)  g_free(maxRMS2);
  if (crossBL1) g_free(crossBL1);
  if (crossBL2) g_free(crossBL2);
  if (corChan)  g_free(corChan);
  if (corIF)    g_free(corIF);
  if (corStok)  g_free(corStok);
  if (BLAnt1)   g_free(BLAnt1);
  if (BLAnt2)   g_free(BLAnt2);
  if (err->error)  Obit_traceback_msg (err, routine, inUV->name);

  /* Give report */
  Obit_log_error(err, OBIT_InfoErr, "%s: flag %8.1lf of %8.1lf vis/interval= %5.1lf percent",
		 routine, (odouble)countBad, (odouble)countAll, 
		 100.0*(odouble)countBad/((odouble)countAll));

  return;
} /* end  ObitUVEditTDRMSAvgVec */

/**
 * Frequency domain editing of visibility data.  
 * Editing is done independently for each visibility measure.  
 * First clipping is done on correlator and Vpol amplitudes.  
 * Following this, an average and RMS is determined for each channel 
 * in each timeAvg period and a spectral baseline is established
 * for the average values, either using a  median window filter (FDwidMW>0) 
 * or a linear baseline fit (FDwidMW<=0) to specified channels.  
 * Channels with excessive RMSes or residual amplitudes are flagged.
 * Flagging is done by entering the offending data in FG table flagTab
 * on outUV.
 * Routine adopted from the AIPSish FDEDIT.FOR/FDEDIT 
 * \param inUV     Name of input uvdata object. 
 *                 Any prior selection and editing is applied.
 * Control parameters on info member of inUV:
 * \li "flagTab" OBIT_int    (1,1,1) FG table version number [ def. 1]
 * \li "timeAvg" OBIT_float  (1,1,1) Time interval over which to average 
 *               data to be flagged (min) [def = 1 min.]
 *               NB: this should be at least 2 integrations.
 * \li "FDmaxAmp"  OBIT_float (1,1,1) Maximum average amplitude allowed in the spectrum 
 *               before fitting.  Any channel exceeding this is flagged in 
 *               advance. default -> infinite 
 * \li "FDmaxV"  OBIT_float (1,1,1) Maximum average amplitude allowed in V polarization; 
 *               any channel exceeding this is flagged in advance of the 
 *               baseline fitting or median filtering, 
 *               Calculates V from difference in amplitudes.   
 *               default -> infinite 
 * \li "FDwidMW"  OBIT_int (1,1,1) If > 0 the width of the median window in channels. 
 *               An odd number (5) is recommended,  default or 0 -> linear baseline
 * \li "FDmaxRMS" OBIT_float (2,1,1) Flag all channels having RMS 
 *               values > maxRMS[0] of the channel median sigma.[default => 6.]
 *               plus maxRMS[1] (default 0.1) of the channel average in quadrature
 * \li "FDmaxRes" OBIT_float (1,1,1) Max. residual flux in sigma allowed for 
 *               channels outside the baseline fitting regions.  
 *               default => 6.
 * \li "FDmaxResBL"  OBIT_float (1,1,1) Max. residual flux in sigma allowed for 
 *               channels within the baseline fitting regions. 
 *               Default = FDmaxRes
 * \li "FDbaseSel"  OBIT_int (4,*,1) Channel selection to define spectral baseline 
 *             Used only for linear baseline fitting.
 *             Select groups of channels/IF(s) to fit as sets 
 *             of (Start,end,inc,IF), i.e., FDbaseSel = 6,37,1,0, 
 *             92,123,1,0 for two regions applying to all IFs.  
 *             Channel and IF numbers 1 -rel
 *             The first group for which the end channel == 0 terminates the list
 *             Channel increments defaults to 1
 *             If the IF==0 then the group applies to all IF.
 *             Default is channels 2 => nchan-1 all IFs
 * \param outUV    UV data onto which the FG table is to be attached.
 *                 May be the same as inUV.
 * \param err      Error stack, returns if not empty.
 */
void ObitUVEditFD (ObitUV* inUV, ObitUV* outUV, ObitErr* err) 
{
  ObitIOCode iretCode, oretCode;
  ObitTableFG *outFlag=NULL;
  ObitTableFGRow *row=NULL;
  gboolean doCalSelect;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  ObitUVDesc *inDesc;
  olong flagTab, iFGRow;
  ofloat timeAvg, maxAmp, maxV, maxResBL, maxRes, maxRMS[2];
  olong  ncorr, numChan, numPol, numIF, numAnt, numBL, widMW, *chanSel;
  olong js, jf, jif, jbl, jndx, indx, kndx, jj,BIF, BChan;
  olong i, j, k, kk, itemp, lastSourceID, curSourceID, lastSubA;
  olong ant1, ant2, blindx, lastFQID=-1;
  ollong  countAll, countBad, *count=NULL;
  olong nVisPIO, ver, firstVis, startVis, suba;
  ofloat startTime, endTime, curTime, amp2, *Buffer;
  ofloat lastTime=-1.0, cbase;
  olong *corChan=NULL, *corIF=NULL, *corStok=NULL, *corV=NULL;
  gboolean *chanMask=NULL, done, gotOne, doMW;
  /* Accumulators per spectral channel/IF/poln/baseline */
  ofloat *sumA=NULL, *sumA2=NULL, *sigma=NULL;
  olong *blLookup=NULL, *BLAnt1=NULL, *BLAnt2=NULL;
  olong defSel[] = {2,-10,1,0, 0,0,0,0};
  ofloat sec;
  gchar *tname, reason[25];
  struct tm *lp;
  time_t clock;
  gchar *routine = "ObitUVEditFD";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));

  /* Fully instantiate UV data if needed */
  ObitUVFullInstantiate (inUV, TRUE, err);
  if (outUV) ObitUVFullInstantiate (outUV, TRUE, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);

  /* Get control parameters */
  flagTab = 1;
  ObitInfoListGetTest(inUV->info, "flagTab", &type, dim, &flagTab);
  /* Time interval */
  timeAvg = 1.0;  /* default 1 min */
  ObitInfoListGetTest(inUV->info, "timeAvg", &type, dim, &timeAvg);
  if (timeAvg<=(1.0/60.0)) timeAvg = 1.0;
  timeAvg /= 1440.0;  /* convert to days */
  /* Maximum amplitude */
  maxAmp = 1.0e20;
  ObitInfoListGetTest(inUV->info, "FDmaxAmp", &type, dim, &maxAmp);
  if (maxAmp<=0.0)  maxAmp = 1.0e20;
  /* Maximum V pol */
  maxV = 1.0e20;  
  ObitInfoListGetTest(inUV->info, "FDmaxV", &type, dim, &maxV);
  if (maxV<=0.0)  maxV = 1.0e20; 
  /* Median window size */
  widMW = 0;  
  ObitInfoListGetTest(inUV->info, "FDwidMW", &type, dim, &widMW);
  /* Maximum average residual */
  maxRes = 6.0;  
  ObitInfoListGetTest(inUV->info, "FDmaxRes", &type, dim, &maxRes);
  if (maxRes<=0.0) maxRes = 6.0; 
  /* Maximum  average residual in baseline fitting area */
  maxResBL = maxRes;  
  ObitInfoListGetTest(inUV->info, "FDmaxResBL", &type, dim, &maxResBL);
  if (maxResBL<=0.0) maxResBL = 6.0; 
  /* RMS clipping parameters */
  maxRMS[0] = 6.0;
  maxRMS[1] = 0.1;
  ObitInfoListGetTest(inUV->info, "FDmaxRMS", &type, dim, maxRMS);
  if ((maxRMS[0]<=0.0) && (maxRMS[1]<=0.0)) {
    maxRMS[0] = 6.0;
    maxRMS[1] = 0.1;
  }
  /* Channel averaging parameters pointer */
  if (!ObitInfoListGetP(inUV->info, "FDbaseSel", &type, dim, (gpointer*)&chanSel)) {
    chanSel = defSel;  /* Use default = channels 2 => n-1 */
  }
  /* FDbaseSel all zero? => default */
  if ((chanSel[0]<=0) && (chanSel[1]<=0) && (chanSel[2]<=0) && (chanSel[3]<=0)) {
    chanSel = defSel;  /* Use default = channels 2 => n-1 */
  }

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

  /* Selection of input? */
  doCalSelect = TRUE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* Open input */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, inUV->name);
  inDesc  = inUV->myDesc;  /* Get descriptor */

  /* Create output FG table */
  tname = g_strconcat ("Flag Table for: ", outUV->name, NULL);
  ver = flagTab;
  outFlag = newObitTableFGValue(tname, (ObitData*)outUV, &ver, OBIT_IO_ReadWrite, 
				err);
  g_free (tname);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
 
  /* Allocate arrays */
  suba    = 1;
  ncorr   = inUV->myDesc->ncorr;
  numAnt  = inUV->myDesc->numAnt[suba-1];/* actually highest antenna number */
  numBL   = (numAnt*(numAnt-1))/2;
  if (inUV->myDesc->jlocs>=0) numPol  = inUV->myDesc->inaxes[inUV->myDesc->jlocs];
  else numPol = 1;
  if (inUV->myDesc->jlocif>=0) numIF  = inUV->myDesc->inaxes[inUV->myDesc->jlocif];
  else numIF = 1;
  if (inUV->myDesc->jlocf>=0) numChan = inUV->myDesc->inaxes[inUV->myDesc->jlocf];
  else numChan = 1;

  /* Make sure it all adds up */
  Obit_return_if_fail ((ncorr==numPol*numIF*numChan), err, 
		       "%s More correlations than Stokes*chan*IF", routine);  

  /* Default baseline fitting? */
  if (chanSel == defSel) defSel[1] = numChan-1;
 
  /* index in order */
  sumA     = g_malloc0 (numChan*numPol*numIF*numBL*sizeof(ofloat));
  sumA2    = g_malloc0 (numChan*numPol*numIF*numBL*sizeof(ofloat));
  sigma    = g_malloc0 (numChan*numPol*numIF*numBL*sizeof(ofloat));
  count    = g_malloc0 (numChan*numPol*numIF*numBL*sizeof(ollong));
  corChan  = g_malloc0 (ncorr * sizeof(olong));     /* Correlator Channel */
  corIF    = g_malloc0 (ncorr * sizeof(olong));     /* Correlator IF */
  corStok  = g_malloc0 (ncorr * sizeof(olong));     /* Correlator Stokes */
  corV     = g_malloc0 (ncorr * sizeof(olong));     /* Stokes L matching R for V */
  chanMask = g_malloc0 (ncorr * sizeof(gboolean)); /* True if channel in baseline region */

  /* Baseline tables */
  blLookup = g_malloc0 (numAnt*sizeof(olong));     /* baseline lookup table */
  BLAnt1   = g_malloc0 (numBL * sizeof(olong));    /* First antenna of baseline */
  BLAnt2   = g_malloc0 (numBL * sizeof(olong));    /* Second antenna of baseline */
  blLookup[0] = 0;
  k = 0;
  for (j=2; j<=numAnt; j++) {BLAnt1[k]=1; BLAnt2[k]=j; k++;}
  for (i=1; i<numAnt; i++) {
    blLookup[i] = blLookup[i-1] + numAnt-i;
    for (j=i+2; j<=numAnt; j++) {BLAnt1[k]=i+1; BLAnt2[k]=j; k++;}
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
  countBad     = 0;

  Buffer = inUV->buffer;

  /* Digest visibility info */
  digestCorrFD (inDesc, corChan, corIF, corStok, corV);

  /* Get mask of channels ordered freq channel, IF, Poln */
  EditFDChanMask (numChan, numIF, numPol, chanSel, chanMask);

  /* Open output table */
  oretCode = ObitTableFGOpen (outFlag, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  
  /* Create Row */
  row = newObitTableFGRow (outFlag);
  
  /* Attach row to output buffer */
  ObitTableFGSetRow (outFlag, row, err);
  if (err->error) goto cleanup;
  
  /* If there are entries in the table, mark it unsorted */
  if (outFlag->myDesc->nrow>0) 
    {outFlag->myDesc->sort[0]=0; outFlag->myDesc->sort[1]=0;}
  
  /* Initialize solution row */
  row->SourID  = 0; 
  row->SubA    = 0; 
  row->freqID  = 0; 
  row->ants[0] = 0; 
  row->ants[1] = 0; 
  row->TimeRange[0] = -1.0e20; 
  row->TimeRange[1] =  1.0e20; 
  row->ifs[0]    = BIF; 
  row->ifs[1]    = 0; 
  row->chans[0]  = BChan; 
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
    
    /* we're in business - loop through data - one vis per read */
    while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
      if ((!gotOne) || (inUV->myDesc->numVisBuff<=0)) { /* need to read new record? */
	if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
	else iretCode = ObitUVRead (inUV, inUV->buffer, err);
	if (err->error) goto cleanup;
	/*if (iretCode!=OBIT_IO_OK) break;*/
	firstVis = inDesc->firstVis;
      }
      
      /* Are we there yet??? */
      done = (inDesc->firstVis >= inDesc->nvis) || 
	(iretCode==OBIT_IO_EOF);
      if (done && (startTime>0.0)) goto process;

      /* Make sure valid data found */
      if (inUV->myDesc->numVisBuff<=0) continue;
      
      gotOne = FALSE;
      /* Time */
      curTime = Buffer[inDesc->iloct];
      if (inDesc->ilocsu>=0) curSourceID = Buffer[inDesc->ilocsu];
      if (startTime < -1000.0) {  /* Set time window etc. if needed */
	startTime = curTime;
	lastTime = curTime;
	endTime   = startTime +  timeAvg;
	startVis  = firstVis;
	lastSourceID = curSourceID;
      }

      /* Still in current interval/source? */
      if ((curTime<endTime) && (curSourceID == lastSourceID) && 
	  (inDesc->firstVis<=inDesc->nvis) && (iretCode==OBIT_IO_OK)) {
	/* accumulate statistics */
	cbase = Buffer[inUV->myDesc->ilocb]; /* Baseline */
	ant1 = (cbase / 256.0) + 0.001;
	ant2 = (cbase - ant1 * 256) + 0.001;
	lastSubA = (olong)(100.0 * (cbase -  ant1 * 256 - ant2) + 1.5);
	/* Baseline index this assumes a1<a2 always - ignore auto correlations */
	if (ant1!=ant2) {
	  blindx =  blLookup[ant1-1] + ant2-ant1-1;
	  if (inDesc->ilocfq>=0) lastFQID = (olong)(Buffer[inDesc->ilocfq]+0.5);
	  else lastFQID = 0;
	  lastTime = curTime;
	  
	  /* Accumulate
	     count =  count per correlation, baseline
	     sumA  = sum of amplitudes for correlation, baseline, then residual
	     sumA  = sum of amplitudes**2 for correlation, baseline, then RMS */
	  indx = inDesc->nrparm; /* offset of start of vis data */
	  for (i=0; i<ncorr; i++) { /* loop 120 */
	    if (Buffer[indx+2] > 0.0) {
	      /* What kind of data is this? */
	      js  = abs(corStok[i])-1;
	      jf  = corChan[i]-1;
	      jif = corIF[i]-1;
	      /* Reorder to freq, IF, poln, baseline */
	      jndx = blindx*ncorr + js*numChan*numIF + jif*numChan + jf;
	      amp2 = Buffer[indx]*Buffer[indx] + Buffer[indx+1]*Buffer[indx+1];
	      count[jndx]++;
	      sumA[jndx]  += sqrt(amp2);
	      sumA2[jndx] += amp2;
	    } 
	    indx += 3;
	  } /* end loop  L120 over correlations */;
	} /* end only cross correlations */
      } else {  /* process interval */
	
      process:
	/* Now have the next record in the IO Buffer */
	if (iretCode==OBIT_IO_OK) gotOne = TRUE;

	/* Average, get RMSes  */
	EditFDAvg(numChan, numIF, numPol, numBL, count, sumA, sumA2, 
		  maxAmp, maxV, corV);
	    
	/* Fit baselines to averages, subtract */
	EditFDBaseline(numChan, numIF, numPol, numBL, count, sumA, sumA2, sigma,
		       widMW, chanMask, err);
	if (err->error) goto cleanup;
	    
	/* Do editing on residuals, RMSes */
	doMW = widMW>0;
	EditFDEdit(numChan, numIF, numPol, numBL, count, sumA, sumA2, sigma,
		   maxRMS, maxRes, maxResBL, chanMask, doMW);

	/* Init Flagging table entry */
	row->SourID  = lastSourceID; 
	row->SubA    = lastSubA; 
	row->freqID  = lastFQID; 
	row->TimeRange[0] = startTime;
	row->TimeRange[1] = lastTime;
	
	/* Loop over baselines/poln/IF/freq, flagging bad */
	indx = 0;
	for (jbl=0; jbl<numBL; jbl++) {
	  for (js=0; js<numPol; js++) {
	    for (jif=0; jif<numIF; jif++) {
	      for (jf=0; jf<numChan; jf++) {

		countAll++;  /* How many examined */

		/* This one flagged? */
		if (sumA[indx]<-9900.0) {
		  countBad++;  /* Count flagged interval/correlator */
		  
		  /* Check for higher number spectral channels in a contigious 
		     block of bad channels and include in same flag */
		  jj = jf;
		  kndx = indx+1;
		  for (kk = jf+1; kk<numChan; kk++) {
		    /* Higher channel number and to be flagged? */
		    if (sumA[kndx]<-9900.0) { /* This one flagged? */
		      jj = kk;          /* This correlator to be included in flag */
		      countBad++;       /* Increment bad count */
		      sumA[kndx] = 0.0; /* don't consider again */
		    } else break;       /* stop looking at higher channel numbers */
		    kndx++;
		  } /* end loop searching for higher bad channels */
		  row->ants[0]   = BLAnt1[jbl]; 
		  row->ants[1]   = BLAnt2[jbl]; 
		  row->ifs[0]    = BIF + jif; 
		  row->ifs[1]    = BIF + jif; 
		  row->chans[0]  = BChan + jf; 
		  row->chans[1]  = BChan + jj; 
		  row->pFlags[0]=row->pFlags[1]=row->pFlags[2]=row->pFlags[3]=0; 
		  /* bit flag implementation kinda screwy */
		  row->pFlags[0] |= 1<<(js);
		  
		  /* Write */
		  iFGRow = -1;
		  oretCode = ObitTableFGWriteRow (outFlag, iFGRow, row, err);
		  if (err->error) goto cleanup;
		} /* end if flagged */
		  
		indx++;
	      } /* end loop over Channel */
	    } /* end loop over IF */
	  } /* end loop over polarization */
	} /* end loop over baseline */

	/* Are we there yet??? */
	done = (inDesc->firstVis >= inDesc->nvis) || 
	  (iretCode==OBIT_IO_EOF);

	/* Reinitialize things */
	startTime = -1.0e20;
	endTime   =  1.0e20;
	for (i=0; i<ncorr*numBL; i++) {
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
  iretCode = ObitUVClose (inUV, err);
  
  /* Close output table */
  oretCode = ObitTableFGClose (outFlag, err);
  
  /* Reset number of vis per read to original value */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV->info, "nVisPIO", OBIT_long, dim, &nVisPIO);

  /* Cleanup */
  /* Deallocate arrays */
  row     = ObitTableFGRowUnref(row);
  outFlag = ObitTableFGUnref(outFlag);
  if (sumA)     g_free(sumA);
  if (sumA2)    g_free(sumA2);
  if (sigma)    g_free(sigma);
  if (count)    g_free(count);
  if (blLookup) g_free(blLookup);
  if (BLAnt1)   g_free(BLAnt1);
  if (BLAnt2)   g_free(BLAnt2);
  if (corChan)  g_free(corChan);
  if (corIF)    g_free(corIF);
  if (corStok)  g_free(corStok);
  if (corV)     g_free(corV);
  if (chanMask) g_free(chanMask);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* Give report */
  Obit_log_error(err, OBIT_InfoErr, "EditFD: flag %8.1lf of %8.1lf vis/interval= %5.1lf percent",
		 (odouble)countBad, (odouble)countAll, 
		 100.0*(odouble)countBad/((odouble)countAll));

  return;
} /* end  ObitUVEditTD */

/**
 * Stokes editing of UV data, FG table out
 *    All data on a given baseline/correlator are flagged if the 
 * amplitude of the datatype "FlagStok"  exceeds maxAmp[0].  
 * If a fraction of bad baselines on any antenna/channel/IF exceeds 
 * maxBad, then all data to that correlator is flagged.  
 *    The actual clipping level is the lesser of maxAmp[0]
 * and if maxAmp[1]<=0, a value determined from a statistical analysis of 
 * each interval intended to flag the most discrepant 3 percent of the data.
 * Flagging entries are written into FG table flagTab.
 * Results are unpredictable for uncalibrated data.
 * Control parameters on info member of inUV:
 * \li "flagStok" OBIT_string (1,1,1) Stokes value to clip (I, Q, U, V, R, L)
 *                default = "V" or correlator based, i.e. RR, LR, XY, YY...
 * \li "flagTab" OBIT_int    (1,1,1) FG table version number [ def. 1]
 *               NB: this should not also being used to flag the input data!
 * \li "timeAvg" OBIT_float  (1,1,1) Time interval over which to determine
 *               data to be flagged (min) [def = 1 min.]
 * \li "maxAmp"  OBIT_float (2,1,1) Maximum value allowed
 *               [0] = clipping level, [1]==0.0 use statistical level
 * \li "minAmp"    OBIT_float (1,1,1) Minimum amplitude allowed
 * \li "maxBad"  OBIT_float (1,1,1) Fraction of allowed flagged baselines 
 *               to an antenna above which all baselines are flagged.
 *               [default 0.25]
 *
 * \param inUV     Input uv data to edit. 
 * \param outUV    UV data onto which the FG table is to be attached.
 *                 May be the same as inUV.
 * \param err      Error stack, returns if not empty.
 */
void ObitUVEditStokes (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  ObitTableFG *outFlag=NULL;
  ObitTableFGRow *row=NULL;
  olong i, j, k, jj, kk, firstVis, startVis, suba, iFGRow, ver;
  ollong countAll, countBad;
  olong lastSourceID, curSourceID, lastSubA, lastFQID=-1;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean gotOne;
  ObitIOAccess access;
  ObitUVDesc *inDesc;
  ofloat timeAvg, lastTime=-1.0, temp[2], maxAmp, maxBad, cbase, hisicc;
  olong *hissig=NULL, *blLookup=NULL;
  olong *corChan=NULL, *corIF=NULL, *corStok=NULL;
  olong *BLAnt1=NULL, *BLAnt2=NULL, BIF, BChan;
  olong flagTab, indx, jndx, kndx,  kndx2, nVisPIO, itemp, ant1, ant2;
  olong numCell, ncorr, numAnt, numBL, blindx, hicell;
  gboolean doStat, done, isbad, *badAnt=NULL;
  ofloat *acc=NULL, *corCnt=NULL, *corBad=NULL, *hiClip=NULL, *Buffer;
  ofloat startTime, endTime, curTime, avg, avg2, hisinc;
  ofloat minAmp=0.0, minAmp2, mxamp2, sec;
  gchar *tname, Stokes[5], oStokes[5], reason[25];
  struct tm *lp;
  time_t clock;
  gchar *routine = "ObitUVEditStokes";

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
  flagTab = 1;
  ObitInfoListGetTest(inUV->info, "flagTab", &type, dim, &flagTab);
  /* Time interval */
  timeAvg = 1.0;  /* default 1 min */
  ObitInfoListGetTest(inUV->info, "timeAvg", &type, dim, &timeAvg);
  if (timeAvg<=(1.0/60.0)) timeAvg = 1.0;
  timeAvg /= 1440.0;  /* convert to days */
 /* RMS clipping parameters, no default */
  ObitInfoListGet(inUV->info, "maxAmp", &type, dim, temp, err);
  maxAmp = temp[0];
  doStat = ((dim[0]==1) || (temp[1]<=0.0));
  mxamp2 = maxAmp * maxAmp;
  ObitInfoListGetTest(inUV->info, "minAmp", &type, dim,  &minAmp);  
  minAmp2 = minAmp*minAmp;
  /* max. fraction bad baselines */
  maxBad = 0.25;           /* default 0.25 */
  ObitInfoListGetTest(inUV->info, "maxBad", &type, dim,  &maxBad);  
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  /* Clipping Stokes type */
  Stokes[0] = 'V'; Stokes[1] = Stokes[2] = Stokes[3] = ' '; Stokes[4] = 0;
  ObitInfoListGetTest(inUV->info, "flagStok",   &type, dim, Stokes);  

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

  /* Set Stokes to ' ' */
  oStokes[0] = oStokes[1] = oStokes[2] = oStokes[3] = ' '; oStokes[4] = 0;
  ObitInfoListGetTest(inUV->info, "Stokes", &type, dim, oStokes);
  dim[0] = strlen(Stokes); dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV->info, "Stokes", OBIT_string, dim, Stokes);

  /* Selection of input - generally will need to convert */
  access = OBIT_IO_ReadCal;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, inUV->name);
  inDesc  = inUV->myDesc;  /* Get descriptor */

  /* Create output FG table */
  tname = g_strconcat ("Flag Table for: ", outUV->name, NULL);
  ver = flagTab;
  outFlag = newObitTableFGValue(tname, (ObitData*)outUV, &ver, OBIT_IO_ReadWrite, 
				err);
  g_free (tname);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
 
  /* Allocate arrays */
  numCell = 800;  /* number of cells in histogram */
  suba    = 1;
  numAnt  = inUV->myDesc->numAnt[suba-1];/* actually highest antenna number */
  numBL   = (numAnt*(numAnt-1))/2;
  ncorr   = inUV->myDesc->ncorr;
  /* acc index = type + corr * (3) + BL * (3*ncorr)  where BL = 0-rel baseline index */
  acc    = g_malloc (ncorr * 3 * numBL * sizeof(ofloat));
  /* hissig index  = cell + numCell*corr */
  hissig = g_malloc (ncorr * numCell * sizeof(olong));
  hisicc = 0.005;  /* Histogram resolution ??? This seems critical */
  hisinc = hisicc*maxAmp; /* Histogram increment */
  hiClip = g_malloc (ncorr * sizeof(ofloat));
  /* corCntt index  = ant + corr*numAnt , ant = 0-rel ant number */
  corCnt = g_malloc (ncorr * numAnt * sizeof(ofloat));
  /* corBad index  = ant + corr*numAnt , ant = 0-rel ant number */
  corBad = g_malloc (ncorr * numAnt * sizeof(ofloat));
  /* badAnt index  = ant + corr*numAnt , ant = 0-rel ant number */
  badAnt = g_malloc (ncorr * numAnt * sizeof(gboolean));
  corChan  = g_malloc (ncorr * sizeof(olong));    /* Correlator Channel */
  corIF    = g_malloc (ncorr * sizeof(olong));    /* Correlator IF */
  corStok  = g_malloc (ncorr * sizeof(olong));    /* Correlator Stokes */

  /* Baseline tables */
  blLookup = g_malloc0(numAnt*sizeof(olong));
  BLAnt1   = g_malloc (numBL * sizeof(olong));    /* First antenna of baseline */
  BLAnt2   = g_malloc (numBL * sizeof(olong));    /* Second antenna of baseline */
  blLookup[0] = 0;
  k = 0;
  for (j=2; j<=numAnt; j++) {BLAnt1[k]=1; BLAnt2[k]=j; k++;}
  for (i=1; i<numAnt; i++) {
    blLookup[i] = blLookup[i-1] + numAnt-i;
    for (j=i+2; j<=numAnt; j++) {BLAnt1[k]=i+1; BLAnt2[k]=j; k++;}
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
  countBad     = 0;
  for (i=0; i<3*ncorr*numBL; i++) acc[i] = 0.0;

  Buffer = inUV->buffer;

  /* Digest visibility info */
  digestCorr (inDesc, &maxAmp, NULL, NULL, NULL, corChan, corIF, corStok);

  /* Open output table */
  oretCode = ObitTableFGOpen (outFlag, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, outFlag->name);
  
  /* Create Row */
  row = newObitTableFGRow (outFlag);
  
  /* Attach row to output buffer */
  ObitTableFGSetRow (outFlag, row, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* If there are entries in the table, mark it unsorted */
  if (outFlag->myDesc->nrow>0) 
    {outFlag->myDesc->sort[0]=0; outFlag->myDesc->sort[1]=0;}
  
  /* Initialize solution row */
  row->SourID  = 0; 
  row->SubA    = 0; 
  row->freqID  = 0; 
  row->ants[0] = 0; 
  row->ants[1] = 0; 
  row->TimeRange[0] = -1.0e20; 
  row->TimeRange[1] =  1.0e20; 
  row->ifs[0]    = BIF; 
  row->ifs[1]    = 0; 
  row->chans[0]  = BChan; 
  row->chans[1]  = 0; 
  row->pFlags[0] = 1<<0 | 1<<1 | 1<<2 | 1<<3; 
  row->pFlags[1] = 0; 
  row->pFlags[2] = 0; 
  row->pFlags[3] = 0; 
  /* Single correlation flags? */
  if (!strncmp(Stokes, "RR", 2)) row->pFlags[0] = 1<<0;
  if (!strncmp(Stokes, "LL", 2)) row->pFlags[0] = 1<<1;
  if (!strncmp(Stokes, "RL", 2)) row->pFlags[0] = 1<<2;
  if (!strncmp(Stokes, "LR", 2)) row->pFlags[0] = 1<<3;
  if (!strncmp(Stokes, "XX", 2)) row->pFlags[0] = 1<<0;
  if (!strncmp(Stokes, "YY", 2)) row->pFlags[0] = 1<<1;
  if (!strncmp(Stokes, "XY", 2)) row->pFlags[0] = 1<<2;
  if (!strncmp(Stokes, "YX", 2)) row->pFlags[0] = 1<<3;
  /* Cross pol */
  if ((Stokes[0]=='Q') && (inDesc->crval[inDesc->jlocs]< 0.0)) 
    row->pFlags[0] = 1<<2 | 1<<3;
  if ((Stokes[0]=='U') && (inDesc->crval[inDesc->jlocs]< 0.0))
    row->pFlags[0] = 1<<2 | 1<<3;

  /* Reason includes time/date */
  /* Get time since 00:00:00 GMT, Jan. 1, 1970 in seconds. */
  time (&clock);
  /* Convert to  broken-down time. */
  lp = localtime (&clock);
  lp->tm_mon++; /* For some bizzare reason, month is 0-rel */
  if (lp->tm_year<1000)  lp->tm_year += 1900; /* full year */
  sec = (ofloat)lp->tm_sec;
  g_snprintf (reason, 25, "Edit %d/%d/%d %d:%d:%f", 
	      lp->tm_year, lp->tm_mon, lp->tm_mday, 
	      lp->tm_hour, lp->tm_min, sec);
  row->reason    = reason; /* Unique string */
  
  /* Loop over intervals */
  done   = FALSE;
  gotOne = FALSE;
  while (!done) {
    
    /* we're in business - loop through data - one vis per read */
    while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
      if ((!gotOne) || (inUV->myDesc->numVisBuff<=0)) { /* need to read new record? */
	iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
	if (err->error) goto cleanup;
	/*if (iretCode!=OBIT_IO_OK) break;*/
	firstVis = inDesc->firstVis;
      }

      /* Are we there yet??? */
      done = (inDesc->firstVis >= inDesc->nvis) || 
	(iretCode==OBIT_IO_EOF);
      if (done && (startTime>0.0)) goto process;

      /* Make sure valid data found */
      if (inUV->myDesc->numVisBuff<=0) continue;

      gotOne = FALSE;
      /* Time */
      curTime = Buffer[inDesc->iloct];
      if (inDesc->ilocsu>=0) curSourceID = Buffer[inDesc->ilocsu];
      if (startTime < -1000.0) {  /* Set time window etc. if needed */
	startTime = curTime;
	lastTime = curTime;
	endTime   = startTime +  timeAvg;
	startVis  = firstVis;
	lastSourceID = curSourceID;
      }

      /* Still in current interval/source? */
      if ((curTime<endTime) && (curSourceID == lastSourceID) && 
	  (inDesc->firstVis<=inDesc->nvis) && (iretCode==OBIT_IO_OK)) {
	/* accumulate statistics */
	cbase = Buffer[inUV->myDesc->ilocb]; /* Baseline */
	ant1 = (cbase / 256.0) + 0.001;
	ant2 = (cbase - ant1 * 256) + 0.001;
	/* DEBUG 
	   if ((ant1==11) && (ant2==28)) {
	   indx = inDesc->nrparm;
	   fprintf (stderr,"time %f vis %f %f %f %f %f %f \n",
	   curTime*24.0,
	   Buffer[indx], Buffer[indx+1], Buffer[indx+2], 
	   Buffer[indx+3], Buffer[indx+4]);
	   }*/
	lastSubA = (olong)(100.0 * (cbase -  ant1 * 256 - ant2) + 1.5);
	/* Baseline index this assumes a1<a2 always - ignore auto correlations */
	if (ant1!=ant2) {
	  blindx =  blLookup[ant1-1] + ant2-ant1-1;
	  if (inDesc->ilocfq>=0) lastFQID = (olong)(Buffer[inDesc->ilocfq]+0.5);
	  else lastFQID = 0;
	  lastTime = curTime;
	  
	  /* Accumulate
	     (1,*) =  count 
	     (2,*) =  sum real then amplitude
	     (3,*) =  sum imaginary then limit */
	  indx = inDesc->nrparm; /* offset of start of vis data */
	  for (i=0; i<ncorr; i++) { /* loop 120 */
	    if (Buffer[indx+2] > 0.0) {
	      jndx = i*3 + blindx*3*ncorr;
	      acc[jndx]   += 1.0;
	      acc[jndx+1] += Buffer[indx];
	      acc[jndx+2] += Buffer[indx+1];
	    } 
	    indx += 3;
	  } /* end loop  L120 over correlations */;
	} /* end only cross correlations */
      } else {  /* process interval */

      process:
	/* Now have the next record in the IO Buffer */
	if (iretCode==OBIT_IO_OK) gotOne = TRUE;

	/* Get average amplitude/collect histogram */
	for (i=0; i<ncorr * numCell; i++)  hissig[i] = 0;
	for (i=0; i<numBL; i++) { /* loop 170 */
	  for (j=0; j<ncorr; j++) { /* loop 160 */
 	    jndx = j*3 + i*3*ncorr;
	    if (acc[jndx] > 0.0) countAll++;  /* count all possibilities */
	    if (acc[jndx] > 0.5) {
	      /* Average amplitude */
	      avg = sqrt (acc[jndx+1]*acc[jndx+1] + acc[jndx+2]*acc[jndx+2]) / acc[jndx];
	      acc[jndx+1] = avg*avg;
	      /* Histogram */
	      hicell = avg / hisinc;
	      /* Could have wrapped */
	      if (hicell<0) hicell = numCell-1;
	      hicell = MIN (hicell, numCell-1);
	      hissig[hicell+j*numCell]++;
	    } 
	  } /* end loop  L160: */
	} /* end loop  L170: */


	/* Histogram flagging?. */
	if (doStat) {
	  editHist (ncorr, numCell, hisinc, hissig, hiClip);
	} else {
	  /* Set histogram clipping levels. */
	  for (j=0; j<ncorr; j++) { /* loop 180 */
            hiClip[j] = maxAmp * maxAmp;
	  } /* end loop  L180: */;
	}

	/* initialize counters */
	for (i=0; i<ncorr*numAnt; i++) corCnt[i] = corBad[i] = 0;

	/* Find bad baselines. */
	for (i=0; i<numBL; i++) { /* loop 200 */
	  for (j=0; j<ncorr; j++) { /* loop 190 */
 	    jndx = j*3 + i*3*ncorr;  /* acc index */
 	    kndx  = BLAnt1[i]-1 + j*numAnt;      /* first antenna index */
 	    kndx2 = BLAnt2[i]-1 + j*numAnt;      /* second antenna index */
	    if (acc[jndx] > 0.5) {
	      avg2 = acc[jndx+1];
	      /* Is this one bad? */
	      acc[jndx+2] = MIN (mxamp2, hiClip[j]);
	      isbad = (avg2  >  acc[jndx+2]) || (avg2 < minAmp2);

	      /* Ant/Correlator info */
	      corCnt[kndx]++;
	      corCnt[kndx2]++;
	      if (isbad) {
		/* Make sure it is flagged. */
		acc[jndx+1] = 1.0e20;
		acc[jndx+2] = 0.0;
		corBad[kndx]++;
		corBad[kndx2]++;
	      } 
	    } else if (acc[jndx] > 0.0) {
	      /* Flag correlators without enough data. */
	      acc[jndx+1] = 1.0e20;
	      acc[jndx+2] = 0.0;
	    }
	  } /* end loop  L190: */
	} /* end loop  L200: */
	
	/* Check for bad antenna/correlators */
	for (i=0; i<numAnt; i++) {
	  for (j=0; j<ncorr; j++) {
	    kndx = i + j*numAnt;
	    badAnt[kndx] = FALSE;
	    if (corCnt[kndx] > 1.1) {
	      if ((corBad[kndx]/corCnt[kndx])  >  maxBad) {
		/* Kill antenna/correlator... */
		badAnt[kndx] = TRUE;
	      }
	    }
	  } 
	}
	
	/* Init Flagging table entry */
	row->SourID  = lastSourceID; 
	row->SubA    = lastSubA; 
	row->freqID  = lastFQID; 
	row->TimeRange[0] = startTime;
	row->TimeRange[1] = lastTime;
	
	/* Loop over antennas/correlator flagging bad */
	for (i=0; i<numAnt; i++) { /* loop 200 */
	  for (j=0; j<ncorr; j++) { /* loop 190 */
 	    kndx  = i + j*numAnt;      /* antenna index */
	    if (badAnt[kndx]) {
	      /* Check for higher number spectral channels in a contigious 
		 block of bad channels and include in same flag */
	      jj = j;
	      for (kk = j+1; kk<ncorr; kk++) {
		/* Only interested in same IF, poln */
		if ((corIF[j]!=corIF[kk]) || (corStok[j]!=corStok[kk])) continue;
		indx = i + kk*numAnt;
		/* Higher channel number and to be flagged? */
		if ((corChan[kk]>corChan[jj]) && badAnt[indx]) {
		  jj = kk;              /* This correlator to be included in flag */
		  badAnt[indx] = FALSE; /* don't consider again */
		} else if (corChan[kk]>corChan[jj]) { /* Higher channel number and good? */
		  break;  /* stop looking at higher channel numbers */
		}
	      } /* end loop searching for higher bad channels */
	      row->ants[0]  = i+1;   /* 1-rel antenna number */
	      row->ants[1]  = 0; 
	      row->ifs[0]   = BIF + corIF[j] - 1; 
	      row->ifs[1]   = BIF + corIF[j] - 1; 
	      row->chans[0] = BChan + corChan[j] - 1; 
	      row->chans[1] = BChan + corChan[jj] - 1; 
	      
	      /* Write */
	      iFGRow = -1;
	      oretCode = ObitTableFGWriteRow (outFlag, iFGRow, row, err);
	      if (err->error) Obit_traceback_msg (err, routine, outFlag->name);
	    } /* end bad correlator section */
	  } /* end loop flagging correlators */
	} /* end loop over antennas */
	
	/* Loop over baselines/correlator flagging bad */
	for (i=0; i<numBL; i++) {
	  for (j=0; j<ncorr; j++) {
	    jndx = j*3 + i*3*ncorr;  /* acc index */
 	    kndx  = BLAnt1[i]-1 + j*numAnt;      /* first antenna index */
 	    kndx2 = BLAnt2[i]-1 + j*numAnt;      /* second antenna index */
	    isbad = badAnt[kndx] || badAnt[kndx2];
	    /* Count flagged interval/correlator */
	    if ((acc[jndx]>0.0) && (isbad || (acc[jndx+1] > acc[jndx+2])))
	      countBad++;
	    if (!isbad) {  /* Don't flag if ant/correlator already flagged */
	      if ((acc[jndx]>0.0) && (acc[jndx+1] > acc[jndx+2])) {
		/* Check for higher number spectral channels in a contigious 
		   block of bad channels and include in same flag */
		jj = j;
		for (kk = j+1; kk<ncorr; kk++) {
		  /* Only interested in same IF, poln */
		  if ((corIF[j]!=corIF[kk]) || (corStok[j]!=corStok[kk])) continue;
		  kndx = kk*3 + i*3*ncorr;
		  /* Higher channel number and to be flagged? */
		  if ((corChan[kk]>corChan[jj]) && 
		      ((acc[kndx]>0.0) && (acc[kndx+1]>acc[kndx+2]))) {
		    jj = kk;         /* This correlator to be included in flag */
		   /* Already counted? countBad++;      Increment bad count */
		    acc[kndx] = 0.0; /* don't consider again */
		  } else if (corChan[kk]>corChan[jj]) { /* Higher channel number and good? */
		    break;  /* stop looking at higher channel numbers */
		  }
		} /* end loop searching for higher bad channels */
		row->ants[0]   = BLAnt1[i]; 
		row->ants[1]   = BLAnt2[i]; 
		row->ifs[0]    = BIF + corIF[j] - 1; 
		row->ifs[1]    = BIF + corIF[j] - 1; 
		row->chans[0]  = BChan + corChan[j]  - 1; 
		row->chans[1]  = BChan + corChan[jj] - 1; 
		
		/* Write */
		iFGRow = -1;
		oretCode = ObitTableFGWriteRow (outFlag, iFGRow, row, err);
		if (err->error) Obit_traceback_msg (err, routine, outFlag->name);
	      }  /* end flag correlator */
	    } /* end correlator not flagged */
	  } /* end loop over correlators */
	} /* end loop over baselines */
	
	/* Are we there yet??? */
	done = (inDesc->firstVis >= inDesc->nvis) || 
	  (iretCode==OBIT_IO_EOF);

	/* Reinitialize things */
	startTime = -1.0e20;
	endTime   =  1.0e20;
	for (i=0; i<3*ncorr*numBL; i++) acc[i] = 0.0;

      } /* end process interval */
      
    } /* end loop processing data */
    if (done) break;
  } /* end loop over intervals */
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) Obit_traceback_msg (err, routine,inUV->name);
  
  /* close uv file */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, inUV->name);
  
  /* Close output table */
  oretCode = ObitTableFGClose (outFlag, err);
  if (err->error) Obit_traceback_msg (err, routine, outFlag->name);
  
  /* Reset number of vis per read to original value */
  dim[0] = 1;
  ObitInfoListAlwaysPut(inUV->info, "nVisPIO", OBIT_long, dim, &nVisPIO);
  /* Reset Stokes */
  dim[0] = strlen(oStokes);
  ObitInfoListAlwaysPut(inUV->info, "Stokes", OBIT_string, dim, oStokes);

  /* Cleanup */
 cleanup:
  /* Deallocate arrays */
  row     = ObitTableFGRowUnref(row);
  outFlag = ObitTableFGUnref(outFlag);
  g_free(acc);
  g_free(hissig);
  g_free(hiClip);
  g_free(corCnt);
  g_free(corBad);
  g_free(badAnt);
  g_free(blLookup);
  g_free(corChan);
  g_free(corIF);
  g_free(corStok);
  g_free(BLAnt1);
  g_free(BLAnt2);

  /* Give report */
  Obit_log_error(err, OBIT_InfoErr, "Edit%c: flag %8.1lf of %8.1lf vis/interval= %5.1lf percent",
		 Stokes[0], (odouble)countBad, (odouble)countAll, 
		 100.0*(odouble)countBad/((odouble)countAll));

  return;
} /* end  ObitUVEditStokes */

/**
 * Clip a uv data set.  Writes edited UV data
 * Control parameters are on the inUV info member:
 * \li "maxAmp" OBIT_float  (1,1,1) Maximum allowed amplitude
 * \li "oper"   OBIT_string (4,1,1) operation type:
 *            "flag" flag data with amplitudes in excess of maxAmp
 *            "clip" clip amplitudes at maxAmp and preserve phase
 *            default is "flag"
 * \param inUV     Input uv data to clip. 
 * \param scratch  True if scratch file desired, will be same type as inUV.
 * \param outUV    If not scratch, then the previously defined output file
 *                 May be the same as inUV.
 *                 May be NULL for scratch only
 *                 If it exists and scratch, it will be Unrefed
 * \param err      Error stack, returns if not empty.
 * \return the clipped ObitUV.
 */
ObitUV* ObitUVEditClip (ObitUV *inUV, gboolean scratch, ObitUV *outUV, 
			ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gboolean doCalSelect;
  gchar *exclude[]={"AIPS CL","AIPS SN","AIPS FG","AIPS CQ","AIPS WX",
		    "AIPS AT","AIPS CT","AIPS OB","AIPS IM","AIPS MC",
		    "AIPS PC","AIPS NX","AIPS TY","AIPS GC","AIPS HI",
		    "AIPS PL",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  olong i, j, indx, firstVis;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ofloat maxAmp, maxAmp2, amp, amp2, temp[2], ratio;
  gboolean same=FALSE, doFlag;
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  gchar oper[5];
  gchar *routine = "ObitUVEditClip";
 
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outUV;
  g_assert (ObitUVIsA(inUV));
  if (!scratch && (outUV==NULL)) {
    Obit_log_error(err, OBIT_Error,"%s Output MUST be defined for non scratch files",
		   routine);
      return outUV;
  }

  /* Fully instantiate UV data if needed */
  ObitUVFullInstantiate (inUV, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outUV);
  if (outUV) ObitUVFullInstantiate (outUV, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, outUV);

  /* Get control parameters */
  /* Maximum amplitude */
  ObitInfoListGet(inUV->info, "maxAmp", &type, dim, temp, err);
  maxAmp = temp[0];
  maxAmp2 = MAX (maxAmp*maxAmp, 1.0e-20);  /* get square */
  oper[0] = 'f'; oper[1] = 'l'; oper[2] = 'a'; oper[3] = 'g'; oper[4] = 0; 
  ObitInfoListGetTest(inUV->info, "oper",   &type, (gint32*)dim,  oper);  
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outUV);
  doFlag = !strncmp (oper, "flag", 4);

 /* Create scratch? */
  if (scratch) {
    if (outUV) outUV = ObitUVUnref(outUV);
    outUV = newObitUVScratch (inUV, err);
  } else if (!same) { /* non scratch output must exist - clone from inUV */
    ObitUVClone (inUV, outUV, err);
  }
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outUV);

  /* Selection of input? */
  doCalSelect = TRUE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* Are input and output the same file? */
  same = ObitUVSame(inUV, outUV, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outUV);

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, inUV->name, outUV);

  /* copy Descriptor */
  if (!same) outUV->myDesc = ObitUVDescCopy(inUV->myDesc, outUV->myDesc, err);

  /* use same data buffer on input and output 
     so don't assign buffer for output */
  if (!same && outUV->buffer) ObitIOFreeBuffer(outUV->buffer);
  outUV->buffer = inUV->buffer;
  outUV->bufferSize = -1;

  /* test open output */
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset output buffer (may be multiply deallocated) */
    outUV->buffer = NULL;
    outUV->bufferSize = 0;
    Obit_traceback_val (err, routine, outUV->name, outUV);
  }

  /* Copy tables before data if different files */
  if (!same) {
    iretCode = ObitUVCopyTables (inUV, outUV, exclude, NULL, err);
    /* If multisource out then copy SU table, multiple sources selected or
       sources deselected suggest MS out */
    if ((inUV->mySel->numberSourcesList>1) || (!inUV->mySel->selectSources))
      iretCode = ObitUVCopyTables (inUV, outUV, NULL, sourceInclude, err);
    if (err->error) {
      outUV->buffer = NULL;
      outUV->bufferSize = 0;
      Obit_traceback_val (err, routine, inUV->name, outUV);
    }
  } /* end of copy tables */

  /* reset to beginning of uv data */
  iretCode = ObitIOSet (inUV->myIO,  inUV->info, err);
  oretCode = ObitIOSet (outUV->myIO, outUV->info, err);
  if (err->error) Obit_traceback_val (err, routine,inUV->name, outUV);

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,inUV->name, outUV);

  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,inUV->name, outUV);

  /* Get descriptors */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;

  /* we're in business, copy, zero data, set weight to 1 */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else iretCode = ObitUVRead (inUV, inUV->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
    if (err->error) Obit_traceback_val (err, routine,inUV->name, outUV);
    firstVis = inUV->myDesc->firstVis;

   /* How many */
    outDesc->numVisBuff = inDesc->numVisBuff;

    /* Clip data */
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
      indx = i*inDesc->lrec + inDesc->nrparm;
      for (j=0; j<inDesc->ncorr; j++) { /* loop over correlations */
	if (inUV->buffer[indx+2] <= 0.0) continue;  /* flagged? */

	/* is value OK */
	amp2 = inUV->buffer[indx]*inUV->buffer[indx] + 
	  inUV->buffer[indx+1]*inUV->buffer[indx+1];
	if (amp2<maxAmp2) continue;

	if (doFlag) {  /* flag excessive values? */
	  inUV->buffer[indx+2] = -fabs(inUV->buffer[indx+2]);
	} /* end if doFlag */ 

	else { /* clip */
	  amp = sqrt (amp2);
	  ratio = maxAmp / MAX (1.0e-10, amp);
	  inUV->buffer[indx]   *= ratio;
	  inUV->buffer[indx+1] *= ratio;
	  inUV->buffer[indx+2] /= ratio;
	}
      } /* end loop over correlations */
    } /* end loop over visibilities */

    /* Write */
    oretCode = ObitUVWrite (outUV, inUV->buffer, err);
    /* suppress vis number update if rewriting the same file */
    if (same) {
      outUV->myDesc->firstVis = firstVis;
      ((ObitUVDesc*)(inUV->myIO->myDesc))->firstVis = firstVis;
    }
  } /* end loop processing data */
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,inUV->name, outUV);
  
  /* unset input buffer (may be multiply deallocated ;'{ ) */
  outUV->buffer = NULL;
  outUV->bufferSize = 0;
  
  /* close files */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, inUV->name, outUV);
  
  oretCode = ObitUVClose (outUV, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, outUV->name, outUV);
  
  return outUV;
} /* end ObitUVEditClip */

/**
 * Clip a uv data set by amplitudes of a given Stokes. Writes edited UV data
 * Data with amplitudes of the selected stokes
 * in excess of maxAmp are flagged.  Optionally all correlations associated
 * may be flagged.  Stokes conversion as needed for test.
 * Control parameters are on the inUV info member:
 * \li "flagStok" OBIT_string (1,1,1) Stokes value to clip (I, Q, U, V, R, L)
 *                default = "I"
 * \li "flagAll"  Obit_bool   (1,1,1) if true, flag all associated correlations
 *                default = True
 * \li "maxAmp"   OBIT_float  (1,1,1) Maximum allowed amplitude
 *
 * \param inUV     Input uv data to clip. 
 * \param scratch  True if scratch file desired, will be same type as inUV.
 * \param outUV    If not scratch, then the previously defined output file
 *                 May be the same as inUV.
 *                 May be NULL for scratch only
 *                 If it exists and scratch, it will be Unrefed
 * \param err      Error stack, returns if not empty.
 * \return the clipped ObitUV.
 */
ObitUV* ObitUVEditClipStokes (ObitUV *inUV, gboolean scratch, ObitUV *outUV, 
			      ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gboolean doCalSelect;
  gchar *exclude[]={"AIPS CL","AIPS SN","AIPS FG","AIPS CQ","AIPS WX",
		    "AIPS AT","AIPS CT","AIPS OB","AIPS IM","AIPS MC",
		    "AIPS PC","AIPS NX","AIPS TY","AIPS GC","AIPS HI",
		    "AIPS PL",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  olong i, indx, cntAll, cntFlagged, firstVis;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ofloat temp[2], maxAmp, maxAmp2, amp2, selFact[2], vis[2], sf1, sf2;
  ofloat wt1, wt2, *Buffer;
  gboolean same, flagAll, flagIt, bothCorr, doConjg;
  olong ichan, iif, istok, nchan, nif, nstok, kstoke0;
  olong incs, incf, incif, jadr[2], ioff, lfoff, ivoff, ip1, ip2;
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  gchar stok[5];
  gchar *routine = "ObitUVEditClipStokes";
 
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outUV;
  g_assert (ObitUVIsA(inUV));
  if (!scratch && (outUV==NULL)) {
    Obit_log_error(err, OBIT_Error,"%s Output MUST be defined for non scratch files",
		   routine);
      return outUV;
  }

  /* Fully instantiate UV data if needed */
  ObitUVFullInstantiate (inUV, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outUV);
  if (outUV) ObitUVFullInstantiate (outUV, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, outUV);

  /* Get control parameters */
  /* Flagging Stokes type */
  stok[0] = 'I'; stok[1] = 0;
  ObitInfoListGetTest(inUV->info, "flagStok",   &type, (gint32*)dim,  stok);  
  flagAll = TRUE;
  ObitInfoListGetTest(inUV->info, "flagAll ",   &type, (gint32*)dim,  &flagAll);  
  /* Maximum amplitude */
  ObitInfoListGet(inUV->info, "maxAmp", &type, dim, temp, err);
  maxAmp = temp[0];
  maxAmp2 = MAX (maxAmp*maxAmp, 1.0e-20);  /* get square */
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outUV);

   /* Are input and output the same file? */
  same = ObitUVSame(inUV, outUV, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outUV);

  /* Create scratch? */
  if (scratch) {
    if (outUV) outUV = ObitUVUnref(outUV);
    outUV = newObitUVScratch (inUV, err);
  } else if (!same) { /* non scratch output must exist - clone from inUV */
    ObitUVClone (inUV, outUV, err);
  }
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outUV);

  /* Selection of input? */
  doCalSelect = TRUE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* Are input and output the same file? */
  same = ObitUVSame(inUV, outUV, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outUV);

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, inUV->name, outUV);

  /* copy Descriptor */
  if (!same) outUV->myDesc = ObitUVDescCopy(inUV->myDesc, outUV->myDesc, err);

  /* use same data buffer on input and output 
     so don't assign buffer for output */
  if (!same && outUV->buffer) {
    ObitIOFreeBuffer(outUV->buffer);
    outUV->buffer = inUV->buffer;
    outUV->bufferSize = -1;
  }

  /* test open output */
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset output buffer (may be multiply deallocated) */
    outUV->buffer = NULL;
    outUV->bufferSize = 0;
    Obit_traceback_val (err, routine, outUV->name, outUV);
  }

  /* Copy tables before data if different files */
  if (!same) {
    iretCode = ObitUVCopyTables (inUV, outUV, exclude, NULL, err);
    /* If multisource out then copy SU table, multiple sources selected or
       sources deselected suggest MS out */
    if ((inUV->mySel->numberSourcesList>1) || (!inUV->mySel->selectSources))
      iretCode = ObitUVCopyTables (inUV, outUV, NULL, sourceInclude, err);
    if (err->error) {
      outUV->buffer = NULL;
      outUV->bufferSize = 0;
      Obit_traceback_val (err, routine, inUV->name, outUV);
    }
  } /* end of copy tables */

  /* reset to beginning of uv data */
  iretCode = ObitIOSet (inUV->myIO,  inUV->info, err);
  oretCode = ObitIOSet (outUV->myIO, outUV->info, err);
  if (err->error) Obit_traceback_val (err, routine,inUV->name, outUV);

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,inUV->name, outUV);

  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,inUV->name, outUV);
  Buffer = inUV->buffer; /* pointer to buffer to work from */

  /* Get descriptors */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;

  /* Set up for parsing data */
  nchan = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocf>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else nif = 1;
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  else nstok = 1;
  /* Get first stokes index */
  if (inDesc->crval[inDesc->jlocs]> 0.0) {
    kstoke0 = inDesc->crval[inDesc->jlocs] + 
      (1.0-inDesc->crpix[inDesc->jlocs]) * inDesc->cdelt[inDesc->jlocs] + 0.5;
  } else {
    kstoke0 = inDesc->crval[inDesc->jlocs] + 
      (1.0-inDesc->crpix[inDesc->jlocs]) * inDesc->cdelt[inDesc->jlocs] - 0.5;
  } 
  /* Data expanded to 3 words per vis */
  incs  = 3 * inDesc->incs  / inDesc->inaxes[0];
  incf  = 3 * inDesc->incf  / inDesc->inaxes[0];
  incif = 3 * inDesc->incif / inDesc->inaxes[0];
  bothCorr = FALSE;
  doConjg  = FALSE;
  
  /* Offsets and stuff by desirec Stokes */
  switch (stok[0]) {
  case ' ':  /* Treat blank as 'I' */
  case 'I':
    if (kstoke0 > 0) {  /* data in Stokes */
      jadr[0] = (1-kstoke0) * incs;
      jadr[1] = jadr[0];
      selFact[0] = 1.0;
      selFact[1] = 0.0;
    } else { /* data correlation based */
      jadr[0] = 0;
      jadr[1] = jadr[0] + incs;
      selFact[0] = 0.5;
      selFact[1] = 0.5;
      /* check if only RR or LL and if so use it. */
      if ((nstok < 2)  ||  (kstoke0 != -1)) jadr[1] = jadr[0];
    }
    break;
  case 'Q':
    if (kstoke0 > 0) {  /* data in Stokes */
      jadr[0] = (2-kstoke0) * incs;
      jadr[1] = jadr[0];
      selFact[0] = 1.0;
      selFact[1] = 0.0;
    } else { /* data correlation based */
      bothCorr = TRUE;    /* Need both correlations for this */
      jadr[0] = (3+kstoke0) * incs;
      jadr[1] = jadr[0] + incs;
      selFact[0] = 0.5;
      selFact[1] = 0.5;
    }
    break;
  case 'U':
    if (kstoke0 > 0) {  /* data in Stokes */
      jadr[0] = (3-kstoke0) * incs;
      jadr[1] = jadr[0];
      selFact[0] = 1.0;
      selFact[1] = 0.0;
    } else { /* data correlation based */
      bothCorr = TRUE;    /* Need both correlations for this */
      doConjg  = TRUE;   /* Need to conjugate result */
      jadr[0] = (3+kstoke0) * incs;
      jadr[1] = jadr[0] + incs;
      selFact[0] = -0.5;
      selFact[1] =  0.5;
    }
    break;
  case 'V':
    if (kstoke0 > 0) {  /* data in Stokes */
      jadr[0] = (4-kstoke0) * incs;
      jadr[1] = jadr[0];
      selFact[0] = 1.0;
      selFact[1] = 0.0;
    } else { /* data correlation based */
      bothCorr = TRUE;    /* Need both correlations for this */
      jadr[0] = (1+kstoke0) * incs;
      jadr[1] = jadr[0] + incs;
      selFact[0] =  0.5;
      selFact[1] = -0.5;
    }
    break;
  case 'R':
    if (kstoke0 > 0) {  /* data in Stokes */
      bothCorr = TRUE;    /* Need both correlations for this */
      jadr[0] = (1-kstoke0) * incs;
      jadr[1] = (4-kstoke0) * incs;
      selFact[0] = 0.5;
      selFact[1] = 0.5;
    } else { /* data correlation based */
      jadr[0] = (1+kstoke0) * incs;
      jadr[1] = jadr[0];
      selFact[0] =  1.0;
      selFact[1] =  0.0;
    }
    break;
  case 'L':
    if (kstoke0 > 0) {  /* data in Stokes */
      bothCorr = TRUE;    /* Need both correlations for this */
      jadr[0] = (1-kstoke0) * incs;
      jadr[1] = (4-kstoke0) * incs;
      selFact[0] = -0.5;
      selFact[1] = -0.5;
    } else { /* data correlation based */
      jadr[0] = (2+kstoke0) * incs;
      jadr[1] = jadr[0];
      selFact[0] =  1.0;
      selFact[1] =  0.0;
    }
    break;
  default:
    Obit_log_error(err, OBIT_Error,"%s Unknown Stokes: %c",
		   routine, stok[0]);
    return outUV;
  }; /* end setup switch 0n Stokes */
  
  /* Abs value of selFact */
  sf1 = fabs(selFact[0]);
  sf2 = fabs(selFact[1]);
  cntAll = 0; cntFlagged = 0;  /* counts */
  
  /* we're in business, check, flag data */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else iretCode = ObitUVRead (inUV, inUV->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
    if (err->error) Obit_traceback_val (err, routine,inUV->name, outUV);
    firstVis = inUV->myDesc->firstVis;

    /* How many in buffer? */
    outDesc->numVisBuff = inDesc->numVisBuff;

    /* Loop over buffer */
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
      indx = i*inDesc->lrec + inDesc->nrparm; /* offset of start of vis data */
      lfoff = -incif - incf;

      /* loop over IF */
      for (iif=0; iif<nif; iif++) {
	lfoff = lfoff + incif;
	ioff  = lfoff;

	/* Loop over frequency channel */
	for (ichan=0; ichan<nchan; ichan++) { /* loop 60 */
	  ioff  = ioff + incf;
	  ivoff = indx + ioff;  /* Buffer offset of real part of first Stokes */

	  /* set input visibility indices of correlations tested */
	  ip1 = ivoff + jadr[0];
	  ip2 = ivoff + jadr[1];

	  /* Need at least one good */
	  if (((sf1*Buffer[ip1+2])<0.0) && ((sf2*Buffer[ip2+2])<0.0)) continue;
	  /* Need both? */
	  if (bothCorr && 
	      (((sf1*Buffer[ip1+2])<0.0) || ((sf2*Buffer[ip2+2])<0.0))) continue;

	  cntAll++;  /* total previously unflagged */

	  /* Do conversion */
	  if (Buffer[ip1+2]>0.0) wt1 = selFact[0];
	  else wt1 = 0.0;
	  if (Buffer[ip2+2]>0.0) wt2 = selFact[1];
	  else wt2 = 0.0;
	  vis[0] = Buffer[ip1+0] * wt1 + Buffer[ip2+0] * wt2;
	  vis[1] = Buffer[ip1+1] * wt1 + Buffer[ip2+1] * wt2;
	  amp2 = (vis[0]*vis[0] +  vis[1]*vis[1]) / ((wt1+wt2)*(wt1+wt2));

	  /* Flag? */
	  flagIt = (amp2 > maxAmp2);
	  if (!flagIt) continue;
	  cntFlagged++;   /* number flagged */

	  /* All stokes or just data in test? */
	  if (flagAll) { /* Flag all Stokes */
	    for (istok=0; istok<nstok; istok++) {
	      Buffer[ivoff+istok*incs+2] = -fabs(Buffer[ivoff+istok*incs+2]);
	    }
	  } else {       /* just test data */
	    Buffer[ip1+2] = - fabs(Buffer[ip1+2]);
	    Buffer[ip2+2] = - fabs(Buffer[ip2+2]);
	  }

	} /* End loop over frequency */
      } /* End loop over IF */     
    } /* End loop over visibilities in buffer */
    
    /* Write */
    oretCode = ObitUVWrite (outUV, inUV->buffer, err);
    /* suppress vis number update if rewriting the same file */
    if (same) {
      outUV->myDesc->firstVis = firstVis;
      ((ObitUVDesc*)(inUV->myIO->myDesc))->firstVis = firstVis;
    }
  } /* End loop processing data */
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,inUV->name, outUV);
  
  /* unset input buffer (may be multiply deallocated ;'{ ) */
  outUV->buffer = NULL;
  outUV->bufferSize = 0;
  
  /* close files */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, inUV->name, outUV);
  
  oretCode = ObitUVClose (outUV, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, outUV->name, outUV);

  /* Give accounting */
  Obit_log_error(err, OBIT_InfoErr,
		 "Flagged %d of %d visibilities with %cPol > %f",
		 cntFlagged, cntAll, stok[0], maxAmp);
  
  return outUV;
} /* end ObitUVEditClipStokes */

/**
 * Flag visibilities by a running mean.
 * A running estimate of the amplitude of each baseline/correlation is determined
 * using a quasi median window (alpha is used to control) filter.
 * Amplitudes in excess of flagSig time a robust RMS in the window will be flagged.
 * Data to be flagged are indicated in FG table flagTab on outUV.
 * Control parameters are on the inUV info member:
 * \li "flagTab"  OBIT_int    (1,1,1) FG table version number [ def. 1]
 *                NB: this should not also being used to flag the input data!
 * \li "flagSig"  OBIT_float (1,1,1) Sigma in window for flagging [def 10]
 * \li "alpha"    OBIT_float (1,1,1) controls averaging [0.5]
 *                0 -> 1 = pure boxcar -> pure MWF (alpha of the 
 *                center data samples are ignored and the rest averaged).
 * \li "timeWind" OBIT_float (1,1,1) Window size (days) over which to determine
 *                medians, RMSes [1 min]
 * \li "timeAvg"  OBIT_float (1,1,1) Previous averaging time in min
 *                if defaulted, determined from data.
 * \li "begIF"    OBIT_int    (1,1,1) First IF in data [ def. 1]
 *                This takes into account any previous selection
 * \li "begChan"  OBIT_int    (1,1,1) First channel in data [ def. 1]
 * \li "allChan"  OBIT_boo    (1,1,1) If true, flag all channels [ def. FALSE]
 *
 * \param inUV    Input uv data to edit. Any prior selection/calibration applied.
 * \param outUV   UV data onto which the FG table is to be attached.
 *                Channels, IS and polarizations will be appropriate for outUV
 *                May be the same as inUV.
 * \param err     Error stack, returns if not empty.
 */
void ObitUVEditMedian (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  ObitTableFG *outFlag=NULL;
  ObitTableFGRow *row=NULL;
  UVMednEditFuncArg** args=NULL;
  ObitThread *myThread=NULL;
  ofloat  flagSig, alpha, timeWind, timeAvg;
  olong flagTab, begIF, begChan; 
  gboolean doCalSelect;
  olong i, j, k, firstVis, startVis, suba, ver;
  ollong checked, countAll, countBad;
  olong nThread;
  olong lastSourceID, curSourceID, lastSubA, lastFQID=-1;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  ObitUVDesc *inDesc;
  ofloat lastTime=-1.0, tInt, cbase, tWind=0.0;
  olong *blLookup=NULL;
  olong *BLAnt1=NULL, *BLAnt2=NULL;
  olong indx, jndx;
  olong numCell, ncorr, numAnt, numBL, blindx, ant1, ant2;
  gboolean gotOne, done, scanDone, scanStartDone, newTime, bufferFull, btemp;
  gboolean allChan = FALSE;
  ofloat *times=NULL, *devs=NULL, *amps=NULL, *Buffer, fblank = ObitMagicF();
  ofloat *work=NULL;
  olong *Chan=NULL, *IFs=NULL, *Stoke=NULL, visNo=-1;
  olong itime=0, ntime, numTime, nextTime, it, itDone, ndata, ndevs;
  ofloat startTime, endTime, curTime=0.0, tTime, scanTime=0.0;
  ofloat sec;
  /*gchar tString[25];*/
  gchar *tname, reason[25];
  struct tm *lp;
  time_t clock;
  gchar *routine = "ObitUVEditMedian";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));

  /* Fully instantiate UV data if needed */
  ObitUVFullInstantiate (inUV, TRUE, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  if (outUV) ObitUVFullInstantiate (outUV, TRUE, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);

  /* Get control parameters */
  flagTab = 1;
  ObitInfoListGetTest(inUV->info, "flagTab", &type, dim, &flagTab);
  begIF = 1;
  ObitInfoListGetTest(inUV->info, "begIF", &type, dim, &begIF);
  begChan = 1;
  ObitInfoListGetTest(inUV->info, "begChan", &type, dim, &begChan);
  allChan = FALSE;
  ObitInfoListGetTest(inUV->info, "allChan", &type, dim, &allChan);
  /* Window Time interval */
  timeWind = 1.0;  /* default 1 min */
  ObitInfoListGetTest(inUV->info, "timeWind", &type, dim, &timeWind);
  if (timeWind<=(1.0/60.0)) timeWind = 1.0;
  timeWind /= 1440.0;  /* convert to days */
  /* alpha */
  alpha = 0.5;
  ObitInfoListGetTest(inUV->info, "alpha", &type, dim, &alpha);
  /* flagging sigma */
  flagSig = 10.0;
  ObitInfoListGetTest(inUV->info, "flagSig", &type, dim,  &flagSig);  
  /* No - flagSig = flagSig*flagSig;  to variance */
  /* Previous averaging time */
  timeAvg = 0.0;
  ObitInfoListGetTest(inUV->info, "timeAvg", &type, dim,  &timeAvg);  
  timeAvg /= 1440.0;

  /* How many integration intervals? */
  tInt = MedianUVInt (inUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  /* Not less that prior averaging if given */
  timeAvg = MAX (tInt, timeAvg);
  numTime = 1 + (olong)(0.5+timeWind/tInt);

  /* How many Threads? */
  myThread = newObitThread();
  nThread = MAX (1, ObitThreadNumProc(myThread));

  /* Create Thread object arrays */
  args = g_malloc0(nThread*sizeof(UVMednEditFuncArg*));
  for (i=0; i<nThread; i++) {
    args[i] = g_malloc0(sizeof(UVMednEditFuncArg));
    args[i]->thread = myThread;
    if (nThread>1) args[i]->ithread = i;
    else           args[i]->ithread = -1;
  }

  /* Look at all data */
  btemp = TRUE;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV->info, "passAll", OBIT_bool, dim, &btemp);
  
  /* Selection of input? */
  doCalSelect = TRUE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* Open input */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, inUV->name);
  inDesc  = inUV->myDesc;  /* Get descriptor */

  /* Create output FG table */
  tname = g_strconcat ("Flag Table for: ", outUV->name, NULL);
  ver = flagTab;
  outFlag = newObitTableFGValue(tname, (ObitData*)outUV, &ver, OBIT_IO_ReadWrite, 
				err);
  g_free (tname);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);

  /* Allocate arrays */
  numCell = 800;  /* number of cells in histogram */
  suba    = 1;
  numAnt  = inUV->myDesc->numAnt[suba-1];/* actually highest antenna number */
  numBL   = (numAnt*(numAnt-1))/2;
  ncorr   = inUV->myDesc->ncorr;
  /* Circular amplitude storage (baseline,correlator,time) */
  ndata   = ncorr * numBL * numTime;
  amps    = g_malloc0 (ndata * sizeof(ofloat));
  /* Time of each time slice */
  times   = g_malloc0 (numTime * sizeof(ofloat));
  /* Maximum deviation in sigma per baseline/correlator */
  ndevs   = ncorr * numBL;
  devs    = g_malloc0 (ndevs * sizeof(ofloat));
  /* Work array for median filtering */
  work    = g_malloc0 (nThread*(numTime+3)*sizeof(ofloat));

  /* Baseline tables */
  blLookup = g_malloc0(numAnt*sizeof(olong));
  BLAnt1   = g_malloc0 (numBL * sizeof(olong));    /* First antenna of baseline */
  BLAnt2   = g_malloc0 (numBL * sizeof(olong));    /* Second antenna of baseline */
  blLookup[0] = 0;
  k = 0;
  for (j=2; j<=numAnt; j++) {BLAnt1[k]=1; BLAnt2[k]=j; k++;}
  for (i=1; i<numAnt; i++) {
    blLookup[i] = blLookup[i-1] + numAnt-i;
    for (j=i+2; j<=numAnt; j++) {BLAnt1[k]=i+1; BLAnt2[k]=j; k++;}
  }

  /* determine channel, IF and Poln ranges for correlations 
   Descriptor on outUV->myIO more accurate */
  Chan  = medianChan (inUV->myDesc, (ObitUVDesc*)outUV->myIO->myDesc, begChan);
  IFs   = medianIFs  (inUV->myDesc, (ObitUVDesc*)outUV->myIO->myDesc, begIF);
  Stoke = medianPoln (inUV->myDesc, (ObitUVDesc*)outUV->myIO->myDesc);

  /* Initialize things */
  startTime    = -1.0e20;
  endTime      =  1.0e20;
  newTime      = FALSE;
  itDone       = 0;
  lastSourceID = -1;
  curSourceID  = 0;
  nextTime     = 0;
  countBad     = 0;
  countAll     = 0;
  checked      = 0;
  ntime        = 1;       /* Number of times in circular buffer */
  scanDone = FALSE;       /* Scan not yet finished */
  scanStartDone = FALSE;  /* Beginning of scan not yet processed */
  /* Single source selected? */
  if (inUV->mySel->sources) curSourceID = inUV->mySel->sources[0];
  lastSubA     = 0;
  for (i=0; i<ndata; i++)   amps[i]  = fblank;  /* Blank accumulator */
  for (i=0; i<ndevs; i++)   devs[i]  = fblank;  /* Blank deviation array */
  for (i=0; i<numTime; i++) times[i] = -1.0e20; /* Blank time array */

  /* Input data I/O buffer */
  Buffer = inUV->buffer;

  /* Open output table */
  oretCode = ObitTableFGOpen (outFlag, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  
  /* Create Row */
  row = newObitTableFGRow (outFlag);
  
  /* Attach row to output buffer */
  ObitTableFGSetRow (outFlag, row, err);
  if (err->error) goto cleanup;
  
  /* If there are entries in the table, mark it unsorted */
  if (outFlag->myDesc->nrow>0) 
    {outFlag->myDesc->sort[0]=0; outFlag->myDesc->sort[1]=0;}
  
  /* Initialize flag row */
  row->SourID  = 0; 
  row->SubA    = 0; 
  row->freqID  = 0; 
  row->ants[0] = 0; 
  row->ants[1] = 0; 
  row->TimeRange[0] = -1.0e20; 
  row->TimeRange[1] =  1.0e20; 
  row->ifs[0]    = 1; 
  row->ifs[1]    = 0; 
  row->chans[0]  = 1; 
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
  g_snprintf (reason, 25, "Median %d/%d/%d %d:%d:%f", 
	      lp->tm_year, lp->tm_mon, lp->tm_mday, 
	      lp->tm_hour, lp->tm_min, sec);
  row->reason    = reason; /* Unique string */
  
  /* Loop over intervals */
  done   = FALSE;
  gotOne = FALSE;  /* Need to read new record */
  while (!done) {
    
    /* we're in business - loop through data - one vis per read */
    while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
      if ((!gotOne) || (inUV->myDesc->numVisBuff<=0)) { /* need to read new record? */
	iretCode = ReadOne (inUV, doCalSelect, &Buffer, &visNo, err);
	if (err->error) goto cleanup;
	/*if (iretCode!=OBIT_IO_OK) break;*/
	firstVis = inDesc->firstVis;
      }
      
      /* Are we there yet??? */
      done = (inDesc->firstVis >= inDesc->nvis) || (iretCode==OBIT_IO_EOF);
      if (done && (startTime>0.0)) goto process;

      /* Make sure valid data found */
      if (inUV->myDesc->numVisBuff<=0) continue;
      
      gotOne = FALSE;
      /* Time */
      curTime = Buffer[inDesc->iloct];
      if (inDesc->ilocsu>=0) curSourceID = Buffer[inDesc->ilocsu];
      /* Scan initialization */
      if (startTime < -1000.0) {  /* Set time window etc. as needed */
	startTime = curTime;
	scanTime  = curTime;
	lastTime  = curTime;
	tWind     = curTime;
	endTime   = startTime + timeWind;
	startVis  = firstVis;
	lastSourceID = curSourceID;
	itime     = 0;
	times[0]  = curTime;
     }

      /* Still in current scan? */
      scanDone = (curTime-lastTime>timeWind)|| /* Gap bigger than window */
	(curSourceID!=lastSourceID)         ||  /* New source */
	(inDesc->firstVis>inDesc->nvis)     ||  /* End of data */
	(iretCode!=OBIT_IO_OK);                 /* Something went wrong in IO */
      if (!scanDone) {
	
	/* Still in current interval? */
	newTime = curTime>tWind;
	
      } else {gotOne = TRUE;  /* Scan done - save this vis for next time */
      } /* End of if scan not done */
      
      /* if buffer full (with timeWind of data ) or scan done
	 - if not done beginning ,do times up to current center
	 - if scan done do all time slices not already done
	 - else do center time
      */
    process:
      /* Buffer just filled with time and need to do start of scan or scan finished or new time */
      bufferFull = ((curTime-scanTime)>=timeWind) && (!scanStartDone);
      if (bufferFull                  || /* Buffer filled  */
	  (newTime && scanStartDone)  || /* Subsequent new time */
	  (iretCode==OBIT_IO_EOF)     || /* Finished file */
	  scanDone) {                    /* finished scan */
	
	/* Init Flagging table entry */
	row->SourID  = lastSourceID; 
	row->SubA    = lastSubA; 
	row->freqID  = lastFQID; 
	
	/* Has the beginning of the scan been done? */
	if ((!scanStartDone) && (!scanDone)) {  /* Beginning of scan */
	  itDone = ntime/2;    /* only do first half */
	  /* Loop over times */
	  for (it=0; it<=itDone; it++) {
	    checked  += MedianDev (amps, it, times, numBL, ncorr, numTime, ntime, alpha, devs, work,
				   nThread, args, err);
	    countBad += MedianFlag (devs, flagSig, numBL, ncorr, allChan, times[it], timeAvg, 
				    BLAnt1, BLAnt2, Chan, IFs, Stoke, outFlag, row, err);
	    if (err->error) goto cleanup;
	  }
	  itDone = it-1;    /* Highest done in current scan */
	  
	  scanStartDone = TRUE;    /* Beginning of scan now processed */
	  
	} else if (scanDone) {  /* All undone intervals in scan */
	  /* Loop over times doing all later than last one done */
	  for (it=0; it<ntime; it++) {
	    if (times[it]>times[itDone]) {
	      checked  += MedianDev (amps, it, times, numBL, ncorr, numTime, ntime, alpha, devs, work,
				     nThread, args, err);
	      countBad += MedianFlag (devs, flagSig, numBL, ncorr, allChan, times[it], timeAvg, 
				      BLAnt1, BLAnt2, Chan, IFs, Stoke, outFlag, row, err);
	      if (err->error) goto cleanup;
	    }
	  }
	  
	  /* Init for next scan */
	  itDone       = 0;
	  startTime    = -1.0e20;
	  nextTime     = 0;
	  ntime        = 1;  /* Number of times in circular buffer */
	  scanStartDone = FALSE;  /* Beginning of next scan not yet processed */
	  for (i=0; i<ndata; i++)   amps[i]  = fblank;  /* Blank accumulator */
	  for (i=0; i<ndevs; i++)   devs[i]  = fblank;  /* Blank deviation array */
	  for (i=0; i<numTime; i++) times[i] = -1.0e20; /* Blank time array */
	  
	} else {  /* just do center time (next undone */
	  itDone++;  /* Highest done in current scan */
	  if (itDone>=ntime) itDone = 0;
	  checked  += MedianDev (amps, itDone, times, numBL, ncorr, numTime, ntime, alpha, devs, work,
				 nThread, args, err);
	  countBad += MedianFlag (devs, flagSig, numBL, ncorr, allChan, times[itDone], timeAvg, 
				  BLAnt1, BLAnt2, Chan, IFs, Stoke, outFlag, row, err);
	  if (err->error) goto cleanup;
	}

	/* Are we there yet??? */
	done = (inDesc->firstVis >= inDesc->nvis) || (iretCode==OBIT_IO_EOF);
	
      } /* end process time slice */
      
      /* Housekeeping for new time slice */
      if (newTime) {
	ntime++;
	ntime = MIN (ntime, numTime);
	tWind = curTime;
	nextTime++;
	newTime = TRUE;
	/* Wrap in circular buffer? */
	if (nextTime>=numTime) nextTime = 0;
	/* Initialize time slot */
	times[nextTime] = curTime;
	/* DEBUG
	T2String (curTime, tString);
	fprintf (stdout,"Time %s\n",tString); */
	/* Blank fill data buffer */
	for (j=0; j<numBL; j++) {
	  for (i=0; i<ncorr; i++) {
	    jndx = numTime * (i + j*ncorr) + nextTime;
	    amps[jndx] = fblank;
	  }
	}
	/* Purge any now expired time slots */
	tTime = curTime-timeWind;
	for (k=0; k<numTime; k++) {
	  if ((times[k]<tTime) && (times[k]>-1.0e10)) {
	    times[k] = -1.0e20;
	    /* Blank fill data buffer */
	    for (j=0; j<numBL; j++) {
	      for (i=0; i<ncorr; i++) {
		jndx = numTime * (i + j*ncorr) + k;
		amps[jndx] = fblank;
	      }
	    }
	  }
	}
	itime = nextTime;
	newTime = FALSE; /* Will deal with this one */
      } /* end new time processing */
      
      /* accumulate data */
      cbase = Buffer[inUV->myDesc->ilocb]; /* Baseline */
      ant1 = (cbase / 256.0) + 0.001;
      ant2 = (cbase - ant1 * 256) + 0.001;
      lastSubA = (olong)(100.0 * (cbase -  ant1 * 256 - ant2) + 1.5);
      /* Baseline index this assumes a1<a2 always - ignore auto correlations */
      if (ant1!=ant2) {
	blindx =  blLookup[ant1-1] + ant2-ant1-1;
	if (inDesc->ilocfq>=0) lastFQID = (olong)(Buffer[inDesc->ilocfq]+0.5);
	else lastFQID = 0;
	lastTime = curTime;

	/* Averaging may cause data not to be exactly in time order 
	   find correct time (within 100 msec) */
	it = itime;
	for (i=0; i<ntime; i++) {
	  if (fabs(curTime-times[i])<1.1574074074074074e-06) {it=i; break;}
	}

	indx = inDesc->nrparm; /* offset of start of vis data */
	for (i=0; i<ncorr; i++) { /* loop 120 */
	  jndx = numTime * (i + blindx*ncorr) + it;
	  if (Buffer[indx+2] > 0.0) {
	    /*amps[jndx] = sqrt (Buffer[indx]*Buffer[indx] + Buffer[indx+1]*Buffer[indx+1]);*/
	    /* Use amp squared */
	    amps[jndx] = (Buffer[indx]*Buffer[indx] + Buffer[indx+1]*Buffer[indx+1]);
	    countAll++;
	  } else {
	    amps[jndx] = fblank;
	  }
	  indx += 3;
	} /* end loop  L120 over correlations */;
 	/* DEBUG 
	if (blindx==0) {
	  fprintf (stdout,"1-2 %d %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f  %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f \n",
		   itime, amps[0], amps[1], amps[2], amps[3], amps[4], amps[5],
		   amps[6], amps[7], amps[8], amps[9], amps[10], amps[11], amps[12]);
	}*/
	
     } /* end only cross correlations */
      /* end accumulate */
      
    } /* end loop read/process data */
    if (done) break;
  } /* end loop over intervals */
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) goto cleanup;
  
  /* Cleanup */
 cleanup:  
  /* Reset passAll */
  btemp = FALSE;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV->info, "passAll", OBIT_bool, dim, &btemp);
  
  /* close uv file */
  iretCode = ObitUVClose (inUV, err);
  
  /* Close output FG table */
  oretCode = ObitTableFGClose (outFlag, err);
  
  /* Deallocate arrays */
  row      = ObitTableFGRowUnref(row);
  outFlag  = ObitTableFGUnref(outFlag);
  if (amps)     g_free(amps);
  if (times)    g_free(times);
  if (devs)     g_free(devs);
  if (work)     g_free(work);
  if (blLookup) g_free(blLookup);
  if (BLAnt1)   g_free(BLAnt1);
  if (BLAnt2)   g_free(BLAnt2);
  if (Chan)     g_free(Chan);
  if (IFs)      g_free(IFs);
  if (Stoke)    g_free(Stoke);
  /* Shut down any threading */
  ObitThreadPoolFree (myThread);
  if (args) {
    for (i=0; i<nThread; i++) if (args[i]) g_free(args[i]);
    g_free(args);
  }
  myThread = ObitThreadUnref(myThread);
  if (err->error)  Obit_traceback_msg (err, routine, inUV->name);

  /* Give report */
  Obit_log_error(err, OBIT_InfoErr, "%s: flag %8.1lf of %8.1lf vis/interval= %6.2lf percent",
		 routine, (odouble)countBad, (odouble)countAll, 
		 100.0*(odouble)countBad/((odouble)countAll));
  if ((100.0*(odouble)checked/((odouble)countAll)) < 0.3) {
    Obit_log_error(err, OBIT_InfoWarn, "%s: Small fraction of data checked (%5.2lf), is timeWind too short?",
		 routine, (odouble)checked/((odouble)countAll));
  }

  return;
} /* end ObitUVEditMedian */

/**
 * Append the contents of the highest numbered AIPS FG table to flagTab
 * And then delete the highest numbered table.
 * Nothing is done if the highest numbered table has a number less than flagTab.
 * \param inUV     Input uv data to with flag tables.
 * \param flagTab  AIPS FG table version number for output, MUST be >= 1
 * \param err      Error stack, returns if not empty.
 */
void ObitUVEditAppendFG (ObitUV *inUV, olong flagTab, ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  ObitTableFG *inFlag=NULL, *outFlag=NULL;
  ObitTableFGRow *row=NULL;
  olong hiTab, iRow, oRow; 
  gchar *routine = "ObitUVEditAppendFG";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));

  Obit_return_if_fail ((flagTab>=1), err, 
		       "%s flagTab MUST be >= 1, not %d", routine, flagTab);  

  /* Fully instantiate UV data if needed */
  ObitUVFullInstantiate (inUV, TRUE, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* Highest table number */
  hiTab = ObitTableListGetHigh (inUV->tableList, "AIPS FG");

  /* If hiTab not greater than flagTab return with warning */
  if (hiTab<=flagTab) {
    Obit_log_error(err, OBIT_InfoWarn, "%s: temp hiTab (%d) must be greater than flagTab (%d)",
		 routine, hiTab, flagTab);
    return;
  }
 
  /* input table */
  inFlag = newObitTableFGValue("InFG", (ObitData*)inUV, &hiTab, OBIT_IO_ReadWrite, 
			       err);
  oretCode = ObitTableFGOpen (inFlag, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;

  /* Anything in input? */
  if (inFlag->myDesc->nrow<=0) goto cleanup;
  
  /* Now output table */
  outFlag = newObitTableFGValue("OutFG", (ObitData*)inUV, &flagTab, OBIT_IO_ReadWrite, 
				err);
  /* Open output table */
  oretCode = ObitTableFGOpen (outFlag, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  
  /* Create Row */
  row = newObitTableFGRow (inFlag);
  
  /* Attach row to output buffer */
  ObitTableFGSetRow (outFlag, row, err);
  if (err->error) goto cleanup;
  
  /* If there are entries in the output table, mark it unsorted */
  if (outFlag->myDesc->nrow>0) 
    {outFlag->myDesc->sort[0]=0; outFlag->myDesc->sort[1]=0;}

  /* Loop over input copying */
  for (iRow=1; iRow<=inFlag->myDesc->nrow; iRow++) {
    iretCode = ObitTableFGReadRow  (inFlag,  iRow, row, err);
    oRow = -1;
    oretCode = ObitTableFGWriteRow (outFlag, oRow, row, err);
    if (err->error) goto cleanup;
  } /* end loop copying */
  
  /* Close tables */
  oretCode = ObitTableFGClose (outFlag, err);
  iretCode = ObitTableFGClose (inFlag, err);

  /* Give report */
  Obit_log_error(err, OBIT_InfoErr, "Copied %d flags from temp table %d to %d",
		 inFlag->myDesc->nrow, hiTab, flagTab);
  /* Cleanup */
 cleanup:
  /* Delete input table */
  iretCode = ObitDataZapTable ((ObitData*)inUV, "AIPS FG", hiTab, err);
  inFlag   = ObitTableFGUnref(inFlag);
  row      = ObitTableFGRowUnref(row);
  outFlag  = ObitTableFGUnref(outFlag);
  if (err->error)  Obit_traceback_msg (err, routine, inUV->name);
} /* end ObitUVEditAppendFG */

/**
 * Routine translated from the AIPSish TDEDIT.FOR/TDCLHI  
 * Determine minimum clipping levels based on histogram of baseline  
 * RMSes.  The clipping level is the median RMS plus 3 times the width  
 * of the RMS distribution or the point which removes 6 percent of the  
 * baselines which ever is greater if normally distributed. 
 * \param ncorr   Number of correlators present 
 * \param numCell Number of cells in each histogram 
 * \param hisinc  Increment of value in histogram 
 * \param hissig  Histogram, index = cell + numCell*corr
 *        On output, the histogram is integrated.
 * \param hiClip  [out] Clipping level in amp**2 units per correlator. 
 */
static void editHist (olong ncorr, olong numCell, ofloat hisinc, olong *hissig, 
		      ofloat* hiClip) 
{
  olong   icorr, count, i, j, ihalf, i6, i16, i56;
  ofloat rhalf, r6, r16, r56, sig, clev;

  if (ncorr == 0) return;
  
  /* Loop over correlator */
  for (icorr= 0; icorr<ncorr; icorr++) { /* loop 400 */
    /* Integrate histogram */
    count = 0;
    for (i=0; i<numCell; i++) { /* loop 50 */
      j = numCell - i - 1;
      count = count + hissig[j+icorr*numCell];
      hissig[j+icorr*numCell] = count;
    } /* end loop  L50:  over cell */
    
    /* Determine levels */
    rhalf = 0.5 * count;
    r6    = 0.06 * count;
    r16   = count * (1.0 / 6.0);
    r56   = count * (5.0 / 6.0);
    ihalf = 1;
    i6    = 1;
    i16   = 1;
    i56   = 2;
    for (i=0; i<numCell; i++) { /* loop 100 */
      if (hissig[i+icorr*numCell] > rhalf) ihalf = i;
      if (hissig[i+icorr*numCell] > r6) i6 = i;
      if (hissig[i+icorr*numCell] > r16) i16 = i;
      if (hissig[i+icorr*numCell] > r56) i56 = i;
    } /* end loop  L100: over cell */

    /* Compute sigma of distribution */
    sig = 0.5 * (i16 - i56) * hisinc;
    
    /* Set clip level */
    clev = MAX ((i6*hisinc), (ihalf*hisinc+3.0*sig));
    
    /* Trap histogram missing */
    if (ihalf >= (numCell-2)) clev = 1000.0;

    /* Return square of flag level */
    hiClip[icorr] = clev * clev;

    /* If the histogram has run off the  high end return a large value */
    if ((i6 >= numCell)  ||  (ihalf >= numCell)) 
      hiClip[icorr] = 10.0;

  } /* end loop  L400: over correlator */
} /* end of routine editHist */ 


/**
 * Determine associated cross polarized correlations, and set maximum 
 * variance for each correlator.  
 * Clip crosspol (maxRMS2) at 2/3 level of parallel.
 * \param inDesc    UV descriptor for data being indexes
 * \param maxRMS    Clipping level, sqrt (maxRMS[0]**2 + (avg_amp * maxRMS[1])**2)
 * \param maxRMS2   [out] maximum variance per correlator,
 *                  in order they occur in the data
 *                  Ignored if NULL.
 * \param crossBL1  [out] if > 0 then index i is for a parallel hand correlation
 *                  with associated cross-hand data correlator crossBL1[i].
 *                  in order they occur in the data
 *                  Ignored if NULL or nstokes <3
 * \param crossBL2  [out] if > 0 then index i is for a parallel hand correlation
 *                  with associated cross-hand data correlator crossBL2[i].
 *                  in order they occur in the data.
 *                  Ignored if NULL or nstokes <3.
 * \param corChan   [out] 1-rel channel number per correlator
 * \param corIF     [out] 1-rel IF number per correlator
 * \param corStok   [out] correlator Stokes -1=RR, -2=LL, -3=RL, -4=LR...
 *                        1=I, 2=Q, 2=U, 4=V
 */
static void digestCorr (ObitUVDesc *inDesc, ofloat *maxRMS, ofloat *maxRMS2, 
			olong *crossBL1, olong *crossBL2, 
			olong *corChan, olong *corIF, olong *corStok) 
{
  olong ichan, iif, istok, nchan, nif, nstok, kstoke0;
  olong incs, incf, incif, jstok, ioff, lfoff, soff;

  /* Set up for parsing data */
  nchan = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocf>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else nif = 1;
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  else nstok = 1;
  /* Get first stokes index */
  if (inDesc->crval[inDesc->jlocs]> 0.0) {
    kstoke0 = inDesc->crval[inDesc->jlocs] + 
      (1.0-inDesc->crpix[inDesc->jlocs]) * inDesc->cdelt[inDesc->jlocs] + 0.5;
  } else {
    kstoke0 = inDesc->crval[inDesc->jlocs] + 
      (1.0-inDesc->crpix[inDesc->jlocs]) * inDesc->cdelt[inDesc->jlocs] - 0.5;
  } 

  /* get increments (one word per correlator) */
  incs  = inDesc->incs  / inDesc->inaxes[0];
  incf  = inDesc->incf  / inDesc->inaxes[0];
  incif = inDesc->incif / inDesc->inaxes[0];

  /* loop over IF */
  for (iif=0; iif<nif; iif++) {
    lfoff = iif * incif;
    ioff  = lfoff;
    
    /* Loop over frequency channel */
    for (ichan=0; ichan<nchan; ichan++) { /* loop 60 */
      soff = ioff;

      /* Loop over polarization */
      for (istok=0; istok<nstok; istok++) {
	corChan[soff]  = ichan+1; 
	corIF[soff]    = iif+1; 

	if (kstoke0>0) { /* Stokes visibility */
	  jstok = kstoke0 + istok;
	  corStok[soff]  = jstok; 
	} else {         /* Correlation */
	  jstok = -kstoke0 + istok;
	  corStok[soff]  = kstoke0 - istok;
	  if ((nstok>2) && (crossBL1!=NULL) && (crossBL2!=NULL)) {
	    /* Is this a correlation with associated cross poln? */
	    if (corStok[soff]==-1) {         /* RR */
	      crossBL1[soff] = soff + 2*incs; 
	      crossBL2[soff] = soff + 3*incs; 
	    } else if (corStok[soff]==-2) {  /* LL */
	      crossBL1[soff] = soff +   incs; 
	      crossBL2[soff] = soff + 2*incs; 
	    } else { /* no */
	      crossBL1[soff] = -1; 
	      crossBL2[soff] = -1; 
	    }
	  } /* end determine cross polarizations */
	}

	/* Set cross polarized clipping at 2/3 of parallel */
	if (maxRMS2!=NULL) {
	  if ((corStok[soff]==-3) || (corStok[soff]==-4)) {
	    maxRMS2[soff] = maxRMS[0] * maxRMS[0] * 0.44444;
	  } else {
	    maxRMS2[soff] = maxRMS[0] * maxRMS[0];
	  }
	}
	
	soff += incs;
      } /* end loop over stokes */
      ioff  += incf;     
    } /* end loop over channel */
  } /* end loop over IF */
  
} /* end digestCorr */

/**
 * Determine associated cross polarized correlations, and set maximum 
 * variance for each correlator.  
 * Clip crosspol (maxRMS2) at 2/3 level of parallel.
 * \param inDesc    UV descriptor for data being indexes
 * \param maxRMSAvg Clipping level RMS/Avg
 * \param maxRMS2   [out] maximum variance per correlator,
 *                  in order they occur in the data
 *                  Ignored if NULL.
 * \param crossBL1  [out] if > 0 then index i is for a parallel hand correlation
 *                  with associated cross-hand data correlator crossBL1[i].
 *                  in order they occur in the data
 *                  Ignored if NULL or nstokes <3
 * \param crossBL2  [out] if > 0 then index i is for a parallel hand correlation
 *                  with associated cross-hand data correlator crossBL2[i].
 *                  in order they occur in the data.
 *                  Ignored if NULL or nstokes <3.
 * \param corChan   [out] 1-rel channel number per correlator
 * \param corIF     [out] 1-rel IF number per correlator
 * \param corStok   [out] correlator Stokes -1=RR, -2=LL, -3=RL, -4=LR...
 *                        1=I, 2=Q, 2=U, 4=V
 */
static void digestCorrTDRMSAvg (ObitUVDesc *inDesc, ofloat maxRMSAvg, ofloat *maxRMS2, 
				olong *crossBL1, olong *crossBL2, 
				olong *corChan, olong *corIF, olong *corStok) 
{
  olong ichan, iif, istok, nchan, nif, nstok, kstoke0;
  olong incs, incf, incif, jstok, ioff, lfoff, soff;

  /* Set up for parsing data */
  nchan = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocf>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else nif = 1;
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  else nstok = 1;
  /* Get first stokes index */
  if (inDesc->crval[inDesc->jlocs]> 0.0) {
    kstoke0 = inDesc->crval[inDesc->jlocs] + 
      (1.0-inDesc->crpix[inDesc->jlocs]) * inDesc->cdelt[inDesc->jlocs] + 0.5;
  } else {
    kstoke0 = inDesc->crval[inDesc->jlocs] + 
      (1.0-inDesc->crpix[inDesc->jlocs]) * inDesc->cdelt[inDesc->jlocs] - 0.5;
  } 

  /* get increments (one word per correlator) */
  incs  = inDesc->incs  / inDesc->inaxes[0];
  incf  = inDesc->incf  / inDesc->inaxes[0];
  incif = inDesc->incif / inDesc->inaxes[0];

  /* loop over IF */
  for (iif=0; iif<nif; iif++) {
    lfoff = iif * incif;
    ioff  = lfoff;
    
    /* Loop over frequency channel */
    for (ichan=0; ichan<nchan; ichan++) { /* loop 60 */
      soff = ioff;

      /* Loop over polarization */
      for (istok=0; istok<nstok; istok++) {
	corChan[soff]  = ichan+1; 
	corIF[soff]    = iif+1; 

	if (kstoke0>0) { /* Stokes visibility */
	  jstok = kstoke0 + istok;
	  corStok[soff]  = jstok; 
	} else {         /* Correlation */
	  jstok = -kstoke0 + istok;
	  corStok[soff]  = kstoke0 - istok;
	  if ((nstok>2) && (crossBL1!=NULL) && (crossBL2!=NULL)) {
	    /* Is this a correlation with associated cross poln? */
	    if (corStok[soff]==-1) {         /* RR */
	      crossBL1[soff] = soff + 2*incs; 
	      crossBL2[soff] = soff + 3*incs; 
	    } else if (corStok[soff]==-2) {  /* LL */
	      crossBL1[soff] = soff +   incs; 
	      crossBL2[soff] = soff + 2*incs; 
	    } else { /* no */
	      crossBL1[soff] = -1; 
	      crossBL2[soff] = -1; 
	    }
	  } /* end determine cross polarizations */
	}

	/* Set cross polarized clipping at 2/3 of parallel */
	if (maxRMS2!=NULL) {
	  if ((corStok[soff]==-3) || (corStok[soff]==-4)) {
	    maxRMS2[soff] = maxRMSAvg * maxRMSAvg * 0.44444;
	  } else {
	    maxRMS2[soff] = maxRMSAvg * maxRMSAvg;
	  }
	}
	
	soff += incs;
      } /* end loop over stokes */
      ioff  += incf;     
    } /* end loop over channel */
  } /* end loop over IF */
  
} /* end digestCorrTDRMSAvg */

/**
 * Determine data type for each correlation fo FD
 * \param inDesc    UV descriptor for data
 * \param corChan   [out] 1-rel channel number per correlator
 * \param corIF     [out] 1-rel IF number per correlator
 * \param corStok   [out] correlator Stokes -1=RR, -2=LL, -3=RL, -4=LR...
 *                        1=I, 2=Q, 2=U, 4=V
 * \param corV      [out] index (0-rel) of Stokes RR correcponding to a Stokes LL
 *                        -1=> none
 */
static void digestCorrFD (ObitUVDesc *inDesc, olong *corChan, olong *corIF, olong *corStok, olong *corV) 
{
  olong ichan, iif, istok, nchan, nif, nstok, kstoke0;
  olong incs, incf, incif, jstok, ioff, lfoff, soff, saveR;

  /* Set up for parsing data */
  nchan = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocf>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else nif = 1;
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  else nstok = 1;
  /* Get first stokes index */
  if (inDesc->crval[inDesc->jlocs]> 0.0) {
    kstoke0 = inDesc->crval[inDesc->jlocs] + 
      (1.0-inDesc->crpix[inDesc->jlocs]) * inDesc->cdelt[inDesc->jlocs] + 0.5;
  } else {
    kstoke0 = inDesc->crval[inDesc->jlocs] + 
      (1.0-inDesc->crpix[inDesc->jlocs]) * inDesc->cdelt[inDesc->jlocs] - 0.5;
  } 

  /* get increments (one word per correlator) */
  incs  = inDesc->incs  / inDesc->inaxes[0];
  incf  = inDesc->incf  / inDesc->inaxes[0];
  incif = inDesc->incif / inDesc->inaxes[0];

  /* loop over IF */
  for (iif=0; iif<nif; iif++) {
    lfoff = iif * incif;
    ioff  = lfoff;
    
    /* Loop over frequency channel */
    for (ichan=0; ichan<nchan; ichan++) { /* loop 60 */
      soff = ioff;
      saveR = -1;

      /* Loop over polarization */
      for (istok=0; istok<nstok; istok++) {
	corChan[soff] = ichan+1; 
	corIF[soff]   = iif+1; 
	corV[soff] = -1;

	if (kstoke0>0) { /* Stokes visibility */
	  jstok = kstoke0 + istok;
	  corStok[soff]  = jstok; 
	} else {         /* Correlation */
	  jstok = -kstoke0 + istok;
	  corStok[soff]  = kstoke0 - istok;
	  if (jstok==1) saveR = soff;
	  if (jstok==2) corV[soff] = saveR;
	}

	soff += incs;
      } /* end loop over stokes */
      ioff  += incf;     
    } /* end loop over channel */
  } /* end loop over IF */
  
} /* end digestCorrFD */

/**
 * Get mask of channels in the linear baseline
 * \param numChan  Number of spectral channels in chanSel, chanMask
 * \param numIF    Number of IFs in chanSel, chanMask
 * \param numPol   Number of polarizations in chanSel, chanMask
 * \param chanSel  Sets of four values (Start,end,inc,IF) giving channels and IFs
 *                 for a block of channels to use in a baseline determination.
 *                 end channel==0 terminates list
 * \param chanMask [out] TRUE if corresponding freq, IF, pol is in baseline
 */
static void EditFDChanMask(olong numChan, olong numIF, olong numPol, olong *chanSel, 
			   gboolean *chanMask)
{
  olong js, jf, jif, indx;
  olong *selTemp;
  gboolean more, match;

  /* Enforce Limits */
  selTemp = chanSel;
  more = selTemp[1]>0;
  while (more) {
    if (selTemp[0]<1)       selTemp[0] = 1;
    if (selTemp[1]>numChan) selTemp[1] = numChan;
    if (selTemp[2]<1)       selTemp[2] = 1;
    if (selTemp[3]<0)       selTemp[3] = 0;
    if (selTemp[3]>numIF)   selTemp[3] = numIF;
    selTemp += 4;
    more = selTemp[1]>0;
  } /* end loop enforcing limits */

  /* Loop over array */
  indx = 0;
  for (js=1; js<=numPol; js++) {
    for (jif=1; jif<=numIF; jif++) {
      for (jf=1; jf<=numChan; jf++) {
	/* Search list to see if this one in it? */
	selTemp = chanSel;
	more = selTemp[1]>0;
	while (more) {
	  /* In IF range? */
	  if ((selTemp[3]<=0) || (selTemp[3]==jif)) {
	    if ((selTemp[0]<=jf) && (selTemp[1]>=jf)) {
	      /* Desired channel? */
	      match = (selTemp[2]<=1) || (((jf-selTemp[0])%(selTemp[2]))==0);
	      chanMask[indx] = match;
	    } /* end channel range */
	  } /* end IF range */
	  selTemp += 4;
	  more = selTemp[1]>0;
	}
	indx++;
      } /* end loop over Channel */
    } /* end loop over IF */
  } /* end loop over polarization */
  
} /* end EditFDChanMask */

/**
 * Average values, get RMSes and do clipping
 * \param numChan Number of spectral channels in count, sumA, sumA2
 * \param numIF   Number of IFs in count, sumA, sumA2
 * \param numPol  Number of polarizations in count, sumA, sumA2
 * \param numBL   Number of baselines in count, sumA, sumA2
 * \param count   Count of entries in each cell of sumA, sumA2
 * \param sumA    [In/Out] Sum of amplitudes in time interval for 
 *                freq, IF, poln, baseline <=0.0 => no data or flagged, 
 *                [Out] average
 *                -9999 => flagged
 * \param sumA2   [In] Sum of amplitude**2 in time interval for 
 *                freq, IF, poln, baseline
 *                [Out] RMS
 * \param maxAmp  Maximum allowable average amplitude
 * \param maxV    Maximum allowable VPol, >1.0e10=> no VPol clipping
 * \param corV    Entry per element of sumA which if >0 is the 0-rel index of the 
 *                corresponding Stokes RR for a Stokes LL measurement, used with maxV
 */
static void EditFDAvg(olong numChan, olong numIF, olong numPol, olong numBL, ollong *count, 
		      ofloat *sumA, ofloat *sumA2, ofloat maxAmp, 
		      ofloat maxV, olong *corV)
{
  olong js, jf, jif, jbl, indx, vindx;
  ofloat vpol;

  /* Loop over array averaging, clip excessive amplitude*/
  indx = 0;
  for (jbl=0; jbl<numBL; jbl++) {
    for (js=0; js<numPol; js++) {
      for (jif=0; jif<numIF; jif++) {
	for (jf=0; jf<numChan; jf++) {
	  if (count[indx]>0) {  /* Any data? */
	    sumA[indx] /= count[indx];
	    if (sumA[indx]>maxAmp) sumA[indx] = -9999.0; /* Too large? */
	    if (count[indx]>1) { /* Enough data for RMS? */
	      sumA2[indx] = sqrt((sumA2[indx]/count[indx] - sumA[indx]*sumA[indx])*
				 ((odouble)count[indx])/(count[indx]-1.0));
	    } else sumA2[indx] = 0.0;
	  } /* end if data */
	  indx++;
	} /* end loop over Channel */
      } /* end loop over IF */
    } /* end loop over polarization */
  } /* end loop over baseline */

  /* Vpol clipping? */
  if (maxV>1.0e10) return;
  /* Loop over array */
  indx = 0;
  for (jbl=0; jbl<numBL; jbl++) {
    vindx = 0;
    for (js=0; js<numPol; js++) {
      for (jif=0; jif<numIF; jif++) {
	for (jf=0; jf<numChan; jf++) {
	  /* Anything to do? */
	  if ((sumA[indx]>0.0) && (corV[vindx]>0) && (sumA[corV[vindx]]>0.0)){ 
	    vpol = fabs(sumA[indx]-sumA[corV[vindx]]);
	    if (vpol>maxV) { /* Flag both */
	      sumA[indx] = sumA[corV[vindx]] = -9999.0;
	    }
	  } /* end if anything to do */
	  vindx++;
	  indx++;
	} /* end loop over Channel */
      } /* end loop over IF */
    } /* end loop over polarization */
  } /* end loop over baseline */

} /* end EditFDAvg */

/**
 * Fit Spectral baseline, either linear or Median Window, subtract
 * Checks than baselines fitted to each IF.
 * \param numChan  Number of spectral channels in count, sumA, RMS
 * \param numIF    Number of IFs in count, sumA, RMS
 * \param numPol   Number of polarizations in count, sumA, RMS
 * \param numBL    Number of baselines in count, sumA, RMS
 * \param count    Count of entries in each cell of avg, RMS
 * \param avg      Sum of amplitudes in time interval for freq, IF, poln, baseline
 *                  <=0.0 => no data or flagged
 * \param RMS      Sum of amplitude**2 in time interval for freq, IF, poln, baseline
 * \param sigma    [out] median alpha sigma for values in avg (MW only)
 * \param widMW    If > 0 the width of the median window in channels.
 *                 If <= 0 use linear baseline and chanMask
 * \param chanMask Mask, True if channel is to be used in baseline fit
 * \param err      Error stack, returns if not empty.
 */
static void EditFDBaseline(olong numChan, olong numIF, olong numPol, olong numBL, 
			   ollong *count, ofloat *avg, ofloat *RMS, ofloat *sigma,
			   olong widMW, gboolean *chanMask, ObitErr* err)
{
  olong js, jf, jif, jbl, jj, jjj, indx, jndx, cnt, half;
  gboolean haveBL, doMW = widMW>0;  /* Median window filter or linear baseline */
  ofloat a, b, aaa, *temp=NULL, *temp2=NULL;
  /*gchar *routine = "EditFDBaseline";*/

  /* error checks */
  if (err->error) return;

  /* Which method? */
  if (doMW) { /* Median window filter */
    temp  = g_malloc0(2*widMW*sizeof(ofloat)); /* Work array for sorting */
    temp2 = g_malloc0(numChan*sizeof(ofloat)); /* Work array for sorting */
    half  = widMW/2;  /* Half width of window */
    /* Loop over array */
    indx = 0;
    for (jbl=0; jbl<numBL; jbl++) {
      for (js=0; js<numPol; js++) {
	for (jif=0; jif<numIF; jif++) {
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
		if ((count[indx+jj]>0) && (temp2[jj]>-9900.0)) {  /* Any data? */
		  temp[cnt++] = temp2[jj];
		}
	      }
	    }
	    /* Get Median and subtract */
	    if ((count[indx]>0) && (avg[indx]>-9900.0)) {
	      aaa = medianVal (temp, 1, cnt);
	      avg[indx] -= aaa;
	      sigma[indx] = MedianSigma (cnt, temp, aaa, 0.2);
	    } else {  /* Datum bad */
	      if (avg[indx]>-9900.0) avg[indx]   = 0.0;
	      sigma[indx] = 100.0;
	    }
	    indx++;
	  } /* end loop over Channel */
	} /* end loop over IF */
      } /* end loop over polarization */
    } /* end loop over baseline */
    g_free(temp);
    g_free(temp2);

  } else { /* Linear baseline fit */
    /* Loop over spectra */
    indx = 0;
    temp = g_malloc0(numChan*sizeof(ofloat)); /* Work array for channel number */
    for (jf=0; jf<numChan; jf++) temp[jf] = (ofloat)jf;
    for (jbl=0; jbl<numBL; jbl++) {
      jndx = 0;
      for (js=0; js<numPol; js++) {
	haveBL = FALSE;
	for (jif=0; jif<numIF; jif++) {

	  /* Fit spectral baseline */
	  EditFDBLFit (temp, &avg[indx], &chanMask[jndx], numChan, &a, &b);
	  jndx += numChan;
	  
	  /* Baseline fitted? */
	  haveBL = a>-900.0;  /* Some baseline data this IF */
	  /* Baseline for this IF? */
	  if (!haveBL) continue;

	  /* Remove baseline */
	  for (jf=0; jf<numChan; jf++) {
	    if (avg[indx]>-9900.0) {
	      avg[indx] -= a + b*temp[jf];
	      indx++;
	    }
	  }
	} /* end loop over IF */
      } /* end loop over polarization */
    } /* end loop over baseline */
   
    g_free(temp);
  } /* end fitting */ 
} /* end  EditFDBaseline */
	    
/**
 * Frequency Domain editing
 * \param numChan  Number of spectral channels in count, avg, RMS
 * \param numIF    Number of IFs in count, avg, RMS
 * \param numPol   Number of polarizations in count, avg, RMS
 * \param numBL    Number of baselines in count, avg, RMS
 * \param count    Count of entries in each cell of avg, RMS
 * \param avg      Sum of amplitudes in time interval for freq, IF, poln, baseline
 *                 <=0.0 => no data or flagged, -9999 => flagged
 * \param RMS      Sum of amplitude**2 in time interval for freq, IF, poln, baseline
 * \param sigma    [out] median alpha sigma for values in avg (MW only)
 * \param maxRMS   Flag all channels having RMS values > maxRMS[0] of the 
 *                 channel median sigma.[default 6.0] plus maxRMS[1] (default 0.1) 
 *                 of the channel average in quadrature
 *                 if maxRMS[0] > 100.0 no RMS flagging
 * \param maxRes   Max. residual flux in sigma allowed for channels outside 
 *                 the baseline fitting regions.
 * \param maxResBL Max. residual flux in sigma allowed for channels within 
 *                 the baseline fitting regions. 
 * \param chanMask Mask, True if channel is to be used in baseline fit
 * \param doMW     True if using median window baselines.
 */
static void EditFDEdit(olong numChan, olong numIF, olong numPol, olong numBL, 
		       ollong *count, ofloat *avg, ofloat *RMS, ofloat *sigma,
		       ofloat *maxRMS, ofloat maxRes, ofloat maxResBL, 
		       gboolean *chanMask, gboolean doMW)
{
  olong js, jf, jif, jbl, indx, jndx, cnt;
  odouble sum1, sum2, meanRMS, residualRMS;
  ofloat *temp=NULL;

  if (maxRMS[0]<100.0)
    temp = g_malloc0(numChan*sizeof(ofloat)); /* Work array for median */

  /* Loop over array */
  indx = 0;
  for (jbl=0; jbl<numBL; jbl++) {
    jndx = 0;
    for (js=0; js<numPol; js++) {
      for (jif=0; jif<numIF; jif++) {
	
	/* Determine average RMS for spectrum */
	if (maxRMS[0]<100.0) {
	  cnt = 0;
	  for (jf=0; jf<numChan; jf++) {
	    if (avg[indx+jf]>-9900.0) temp[cnt++] = RMS[indx+jf];
	  }
	  
	  /* Median rms */
	  if (cnt>0) meanRMS = medianVal(temp, 1, cnt);
	  else meanRMS = 0.0;
	  
	  /* RMS flagging */
	  for (jf=0; jf<numChan; jf++) {
	    if ((avg[indx+jf]>-9900.0) && 
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
	    cnt++; 
	    sum1 += avg[indx+jf];
	    sum2 += avg[indx+jf]*avg[indx+jf];
	  }
	}
	/* Residual RMS */
	if (cnt>0) residualRMS = sqrt (sum2/cnt + (sum1/cnt)*(sum1/cnt));
	else residualRMS = 1.0e20;
	if (!doMW) sigma[indx+jf] = residualRMS;

	/* Residual flagging */
	for (jf=0; jf<numChan; jf++) {
	  if (chanMask[jndx+jf]) { /* In baseline */
	    if (fabs(avg[indx+jf])>maxResBL*sigma[indx+jf]) avg[indx+jf] = -9999.0;
	  } else {  /* Not in baseline */
	    if (fabs(avg[indx+jf])>maxRes*sigma[indx+jf])   avg[indx+jf] = -9999.0;
	  }
	} /* end loop over Channel */
	jndx += numChan;
	indx += numChan;
      } /* end loop over IF */
    } /* end loop over polarization */
  } /* end loop over baseline */

  if (temp) g_free(temp); /* Cleanup */
  
} /* end EditFDEdit */

/**
 * Return median of an array
 * Does sort then returns median
 * \param array array of values, on return will be in ascending order
 * \param incs  increment in data array
 * \param n     dimension of array
 * \return      Median/average, zero if no data
 */
static ofloat medianVal (ofloat *array, olong incs, olong n)
{
  ofloat out=0.0;

  if (n<=0) return out;

  /* Sort to ascending order */
  qsort ((void*)array, n, MIN(2,incs)*sizeof(ofloat), compare_ofloat);

  out = array[n/2];

  return out;
} /* end medianVal */
/**
 * ofloat comparison of two arguments
 * \param arg1 first value to compare
 * \param arg2 second value to compare
 * \return negative if arg1 is less than arg2, zero if equal
 *  and positive if arg1 is greater than arg2.
 */
static int compare_ofloat  (const void* arg1,  const void* arg2)
{
  int out = 0;
  ofloat larg1, larg2;

  larg1 = *(ofloat*)arg1;
  larg2 = *(ofloat*)arg2;
  if (larg1<larg2)      out = -1;
  else if (larg1>larg2) out = 1;
  return out;
} /* end compare_ofloat */

/**
 * Routine to fit Linear slope:  
 * y = a + b * x  where m
 * Routine translated from the AIPSish FDEDIT.FOR/FDFLFT  
 * \param x    "X" values of the points 
 * \param y    "Y" values of the points 
 * \param w     Where true use x,y
 * \param ndata Number of data points 
 * \param a     [out] Zeroth order polynomial coefficient, 
 *              -1000.0 if no data
 * \param b     [out] First order polynomial coefficient 
 */
static void EditFDBLFit (ofloat* x, ofloat* y, gboolean *m, olong ndata, 
			 ofloat *a, ofloat *b) 
{
  odouble sx, sy, sxy, syy, sxx;
  olong  i, count;
  ofloat     delta;

  sx = sy = sxx = sxy = syy = 0.0;
  count = 0;

  /* Make sums */
  for (i=0; i<ndata; i++) { /* loop 100 */
    if (m[i]) {
      count++;
      sy  += y[i];
      sx  += x[i];
      sxx += x[i] * x[i];
      sxy += x[i] * y[i];
      syy += y[i] * y[i];
    } 
  } /* end loop  L100: */

  /* Results */
  (*a) = -1000.0;
  (*b) = 0.0;
  if (count >= 2) {
    delta = count * sxx - sx * sx;
    if (delta > 0) {
      (*a) = ( sxx * sy - sx * sxy ) / delta;
      (*b) = ( sxy * count - sx * sy ) / delta;
    } 
  }
} /* end EditFDBLFit */ 

/**
 * Estimate data integration
 * Measure the minimum change greater than zero in data label time 
 * in first 20 time changes.
 * \param inUV     Input uv data to edit. 
 * \param err      Error stack, returns if  error.
 * \return data integration time in days;

 */
static ofloat MedianUVInt (ObitUV *inUV, ObitErr *err)
{
  ofloat out=1.0e20;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong i, indx, count;
  ObitUVDesc *inDesc;
  ObitIOCode iretCode;
  gboolean doCalSelect;
  ofloat  lastTime, curTime, tInt;
  ObitIOAccess access;
  gchar *routine = "ObitUVEditMedian";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitUVIsA(inUV));

  /* Selection of input? */
  doCalSelect = TRUE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;
  
  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, inUV->name, out);

  inDesc  = inUV->myDesc;  /* Get descriptor */
  lastTime = -1.0e20;      /* Init last time */
  count    = 0;
  tInt     = 1.0e20;       /* Init integration time */

  iretCode = OBIT_IO_OK;
  while (iretCode == OBIT_IO_OK) {
    if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else iretCode = ObitUVRead (inUV, inUV->buffer, err);
    if (err->error)  Obit_traceback_val (err, routine, inUV->name, out);
    for (i=0; i<inUV->myDesc->numVisBuff; i++) { /* loop over buffer */
      indx = i*inDesc->lrec ;
      curTime = inUV->buffer[indx+inDesc->iloct];
      if (curTime>lastTime) {
	if ((curTime-lastTime)>0.0) tInt = MIN (tInt, curTime-lastTime);
	lastTime = curTime;
	count++;
	if (count>20) break;
      }
    } /* end loop over buffer */
  } /* End loop reading data */

  out = tInt;

  /* close uv file */
  iretCode = ObitUVClose (inUV, err);
  if (err->error)  Obit_traceback_val (err, routine, inUV->name, out);

  return out;
} /* end MedianUVInt */

/**
 * Determine deviations in sigma from (alpha) median for a given time 
 * \param amps    Circular amplitude storage (baseline,correlator,time)
 * \param itime   time slice of interest
 * \param times   times of data samples, <-1000 => ignore
 * \param numBL   number of baselines
 * \param numCorr number of correlations
 * \param numTime number of allocated time slices in amps
 * \param ntime   number of actual time slices
 * \param alpha   controls averaging
 *               0 -> 1 = pure boxcar -> pure MWF (alpha of the 
 *               center data samples are ignored and the rest averaged).
 * \param devs   [out] Deviations in sigma per baseline/correlator
 *               fblank => not determined
 * \param work   work array at least the size of ntime per thread
 * \param nThread Number of threads to use
 * \param args    Array of thread arguments
 * \param err     Error stack, returns if  error.
 * \return number of baselines/correlations with valid data
 */
static ollong MedianDev (ofloat *amps, olong itime, ofloat *times,
			 olong numBL, olong numCorr, olong numTime, olong ntime,
			 ofloat alpha, ofloat *devs, 
			 ofloat *work, olong nThread, UVMednEditFuncArg** args,
			 ObitErr *err)
{
  ollong out = 0;
  ObitThreadFunc func=(ObitThreadFunc)ThreadMedianDev;
  olong  i, nBLPerThread, loBL, hiBL;
  gboolean OK;
  gchar *routine = "MedianDev";

  /* error checks */
  if (err->error) return out;

  /* Divide up work */
  nBLPerThread = numBL/nThread;
  loBL = 1;
  hiBL = nBLPerThread;
  hiBL = MIN (hiBL, numBL);
  
  /* Set up thread arguments */
  for (i=0; i<nThread; i++) {
    if (i==(nThread-1)) hiBL = numBL;  /* Make sure do all */
    args[i]->first  = loBL;
    args[i]->last   = hiBL;
    if (nThread>1) args[i]->ithread =  i;
    else           args[i]->ithread = -1;
    /* Other stuff */
    args[i]->amps    = amps;
    args[i]->itime   = itime;  
    args[i]->times   = times;   
    args[i]->numBL   = hiBL-loBL+1;;
    args[i]->numCorr = numCorr; 
    args[i]->numTime = numTime; 
    args[i]->ntime   = ntime;  
    args[i]->alpha   = alpha;   
    args[i]->devs    = devs;   
    args[i]->work    = &work[i*ntime+2];
    args[i]->number  = 0; 
    /* Update which BL range */
    loBL += nBLPerThread;
    hiBL += nBLPerThread;
    hiBL = MIN (hiBL, numBL);
  }
  
  /* Do operation on buffer possibly with threads */
  OK = ObitThreadIterator (args[0]->thread, nThread, func, (gpointer**)args);
  
  /* Check for problems */
  if (!OK) {
    Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
  }
  
  /* Sum counts */
  for (i=0; i<nThread; i++) out += args[i]->number;
  
  return out;
} /* end MedianDev */ 

/**
 * Return a single visibility from a buffer read
 * Modified descriptor values  firstVis and numVisBuff.
 * \param inUV    Input uv data to read. Any prior selection/calibration applied.
 * \param doCalSelect If TRUE apply calibration/selection
 * \param Buffer  [out] pointer to start of vis
 * \param visNo   [in/out] Visibility number (1-rel) in buffer for read
 *                on input the number of the last read
 * \param err     Error stack, returns if  error.
 * \return return code, OBIT_IO_OK => OK
 */
static ObitIOCode ReadOne (ObitUV *inUV, gboolean doCalSelect, ofloat** Buffer, 
			   olong *visNo, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  gchar *routine = "ReadOne";

  /* error checks */
  if (err->error) return retCode;

  /* Data still in buffer? */
  if (((*visNo)>0) && ((*visNo)<=((ObitUVDesc*)inUV->myIO->myDesc)->numVisBuff)) {
    /* Yes */
    *Buffer = inUV->buffer + (*visNo)*inUV->myDesc->lrec;
    inUV->myDesc->firstVis   = ((ObitUVDesc*)inUV->myIO->myDesc)->firstVis + (*visNo);
    inUV->myDesc->numVisBuff = ((ObitUVDesc*)inUV->myIO->myDesc)->numVisBuff - (*visNo);
    (*visNo)++;
    return retCode;
  }
  /* Must read */
  /* Restore state */
  inUV->myDesc->firstVis   = ((ObitUVDesc*)inUV->myIO->myDesc)->firstVis;
  inUV->myDesc->numVisBuff = ((ObitUVDesc*)inUV->myIO->myDesc)->numVisBuff;
  /* Read */
  if (doCalSelect) retCode = ObitUVReadSelect (inUV, inUV->buffer, err);
  else retCode = ObitUVRead (inUV, inUV->buffer, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, retCode);

  /* Set output */
  *Buffer = inUV->buffer;
  *visNo  = 1;
  return retCode;
} /* end ReadOne */
/**
 * Determine alpha median value of a ofloat array
 * \param n       Number of points
 * \param value   Array of values
 * \param alpha   0 -> 1 = pure boxcar -> pure MWF (ALPHA of the 
 *                data samples are discarded and the rest averaged). 
 * \return alpha median value
 */
static ofloat MedianLevel (olong n, ofloat *value, ofloat alpha)
{
  ofloat out=0.0;
  ofloat fblank = ObitMagicF();
  ofloat beta, sum;
  olong i, i1, i2, count;

  if (n<=0) return out;

  /* Sort to ascending order */
  qsort ((void*)value, n, sizeof(ofloat), compare_ofloat);

  out = value[n/2];

  beta = MAX (0.05, MIN (0.95, 1.0-alpha)) / 2.0; /*  Average around median factor */

  /* Average around the center */
  i1 = MAX (0, (n/2)-(olong)(beta*n+0.5));
  i2 = MIN (n, (n/2)+(olong)(beta*n+0.5));

  if (i2>i1) {
    sum = 0.0;
    count = 0;
    for (i=i1; i<i2; i++) {
      if (value[i]!=fblank) {
	sum += value[i];
	count++;
      }
    }
    if (count>0) out = sum / count;
  }
   
  return out;
} /* end MedianLevel */

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
static ofloat MedianSigma (olong n, ofloat *value, ofloat mean, ofloat alpha)
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
} /* end MedianSigma */


/**
 * Write flag table 
 * \param devs    Deviations in sigma per baseline/correlator
 *                fblank => not determined
 * \param flagSig Flagging level in devs
 * \param numBL   Number of baselines
 * \param numCorr Number of correlators
 * \param allChan If TRUE flag all channels
 * \param time    Time of this time interval
 * \param timeInt Width in time to flag
 * \param BLAnt1  First antenna number of baseline
 * \param BLAnt2  Second antenna number of baseline
 * \param Chan    channel range to flag, 2 values per baseline/correlator
 *                in same order as devs
 * \param IFs     IF range to flag, 2 values per baseline/correlator
 *                in same order as devs
 * \param Stoke   Stokes flag, sum of 1=RR, 2=LL, 4=RL, 8=LR
 *                in same order as devs, 
 * \param FGtab   AIPS FG table to write
 * \param FGRow   AIPS FG table row to write
 * \param err     Error stack, returns if error.
 * \return number of baselines/correlations flagged
 */
static ollong MedianFlag (ofloat *devs, ofloat flagSig, 
			  olong numBL, olong numCorr, gboolean allChan,
			  ofloat time, ofloat timeInt,
			  olong *BLAnt1, olong *BLAnt2, 
			  olong *Chan, olong *IFs, olong *Stoke,
			  ObitTableFG *FGtab, ObitTableFGRow *FGrow, 
			  ObitErr *err)
{
  ollong out = 0;
  olong  iBL, icorr, indx, jndx, kndx;
  olong iFGRow;
  ofloat fblank = ObitMagicF();

  gchar *routine = "MedianFlag";

  /* error checks */
  if (err->error) return out;
  if (time<-1000.0) return out;  /* Data for this time? */

  /* Loop over baselines */
  for (iBL=0; iBL<numBL; iBL++) {
    jndx = -2;
    kndx = -1;
    /* Loop over correlator */
    for (icorr=0; icorr<numCorr; icorr++) {
      jndx += 2;
      kndx += 1;
      
      /* Flag it? */
      indx = icorr + iBL*numCorr;
      if ((devs[indx]!=fblank) && (devs[indx]>flagSig)) {
	/* Flag */
	out++;   /* Count flagged data */
	FGrow->TimeRange[0] = time - timeInt*0.5;
	FGrow->TimeRange[1] = time + timeInt*0.5;
	FGrow->ants[0]  = BLAnt1[iBL]; 
	FGrow->ants[1]  = BLAnt2[iBL]; 
	FGrow->ifs[0]   = IFs[jndx];
	FGrow->ifs[1]   = IFs[jndx+1];
	if (allChan) {
	  FGrow->chans[0] = 1;
	  FGrow->chans[1] = 0;
	} else {
	  FGrow->chans[0] = Chan[jndx];
	  FGrow->chans[1] = Chan[jndx+1];
	}
	FGrow->pFlags[0]=FGrow->pFlags[1]=FGrow->pFlags[2]=FGrow->pFlags[3]=0; 
	/* bit flag implementation kinda screwy */
	if (Stoke[kndx]==1) FGrow->pFlags[0] |= 1<<(0);
	else if (Stoke[kndx]==2) FGrow->pFlags[0] |= 1<<(1);
	else if (Stoke[kndx]==4) FGrow->pFlags[0] |= 1<<(2);
	else if (Stoke[kndx]==8) FGrow->pFlags[0] |= 1<<(3);
	/* DEBUG
	fprintf (stdout," %f %d %d %d %f\n",
		 time, BLAnt1[iBL], BLAnt2[iBL], Stoke[kndx], devs[indx]);*/
	/* Write */
	iFGRow = -1;
	ObitTableFGWriteRow (FGtab, iFGRow, FGrow, err);
	if (err->error)  Obit_traceback_val (err, routine, FGtab->name, out);
	devs[indx] = fblank;  /* Don't repeat this one */
      }
    } /* end loop over correlator  */
  } /* end loop over baseline */
  return out;
} /* end MedianFlag */ 

/**
 * Get array of (start, end) channels (1-rel) for correlations in inDesc
 * These are derived assuming the data described by inDesc are
 * described by outDesc.
 * \param inDesc   Descriptor for data being flagged
 * \param outDesc  Descriptor for data before averaging.
 * \param begChan  Actual 1st Channel (1-rel) after prior selection
 * \return array of (start,stop) channels for each correlation in inDesc
 *  Should be g_freed when done.
 */
static olong* medianChan (ObitUVDesc *inDesc, ObitUVDesc *outDesc, olong begChan)
{
  olong *out=NULL;
  olong i, jincs, jincf, jincif, nif;
  olong ichan, iif, ipoln, nchAvg, bch, ech;

  out = g_malloc0(2*inDesc->ncorr*sizeof(olong));  /* Allocate output array */

  /* Increments */
  jincs = inDesc->incs / inDesc->inaxes[0];
  jincf = inDesc->incf / inDesc->inaxes[0];
  if (inDesc->jlocif>0) {
    jincif = inDesc->incif / inDesc->inaxes[0];
    nif    = inDesc->inaxes[inDesc->jlocif];
  }  else {
    jincif = 1;
    nif    = 1;
  }
  
  /* Number of channels averaged */
  nchAvg = outDesc->inaxes[outDesc->jlocf] /  inDesc->inaxes[inDesc->jlocf];
  nchAvg = MAX (1, nchAvg);  
  /* If selection - assume no averaging */
  if (begChan>1) nchAvg = 1;

  /* Loop over channels */
  for (ichan=0; ichan<inDesc->inaxes[inDesc->jlocf]; ichan++) {
    bch = (begChan-1) + ichan*nchAvg;
    ech = bch + nchAvg;
    if (nchAvg==outDesc->inaxes[outDesc->jlocf]) ech = 0;  /* All? */
    if (inDesc->inaxes[outDesc->jlocf]==1) ech = bch+1;
    /* Loop over IF */
    for (iif=0; iif<nif; iif++) {
      /* Loop over poln */
      for (ipoln=0; ipoln<inDesc->inaxes[inDesc->jlocs]; ipoln++) {
	i = 2 * (ichan*jincf + iif*jincif + ipoln*jincs);
	out[i]   = bch + 1;
	out[i+1] = ech;
      }
    }
  }

  return out;
} /* end medianChan */

/**
 * Get array of (start, end) IFs (1-rel) for correlations in inDesc
 * These are derived assuming the data described by inDesc are
 * described by outDesc.
 * \param inDesc   Descriptor for data being flagged
 * \param outDesc  Descriptor for data before averaging.
 * \param begIF    Actual 1st IF (1-rel) after prior selection
 * \return array of (start,stop) IFs for each correlation in inDesc
 *  Should be g_freed when done.
 */
static olong* medianIFs (ObitUVDesc *inDesc, ObitUVDesc *outDesc, olong begIF)
{
  olong *out=NULL;
  olong i, jincs, jincf, jincif, nif;
  olong ichan, iif, ipoln, nifAvg, bif, eif;

  out = g_malloc0(2*inDesc->ncorr*sizeof(olong));  /* Allocate output array */

  /* Increments */
  jincs = inDesc->incs / inDesc->inaxes[0];
  jincf = inDesc->incf / inDesc->inaxes[0];
  if (inDesc->jlocif>0) {
    jincif = inDesc->incif / inDesc->inaxes[0];
    nif    = inDesc->inaxes[inDesc->jlocif];
  }  else {
    jincif = 1;
    nif    = 1;
  }
  
  /* Number of IFs averaged */
  nifAvg = outDesc->inaxes[outDesc->jlocif] /  inDesc->inaxes[inDesc->jlocif];
  nifAvg = MAX (1, nifAvg);  
  /* If selection - assume no averaging */
  if (begIF>1) nifAvg = 1;
 
  /* Loop over IF */
  for (iif=0; iif<nif; iif++) {
    bif = (begIF-1) + iif*nifAvg;
    eif = (begIF-1) + (iif+1)*nifAvg;
    if (nifAvg==outDesc->inaxes[outDesc->jlocif]) eif = 0;  /* All? */
    if (inDesc->inaxes[outDesc->jlocif]==1) eif = bif+1;
    /* Loop over channels */
    for (ichan=0; ichan<inDesc->inaxes[inDesc->jlocf]; ichan++) {
      /* Loop over poln */
      for (ipoln=0; ipoln<inDesc->inaxes[inDesc->jlocs]; ipoln++) {
	i = 2 * (ichan*jincf + iif*jincif + ipoln*jincs);
	out[i]   = bif + 1;
	out[i+1] = eif;
      }
    }
  }
  
  return out;
} /* end  medianIFs */

/**
 * Get array of polarization flags, sum of 1=RR, 2=LL, 4=RL, 8=LR
 * These are derived assuming the data described by inDesc are
 * described by outDesc.
 * \param inDesc   Descriptor for data being flagged
 * \param outDesc  Descriptor for data before averaging.
 * \return array of polarization flags for each correlation in inDesc
 *  Should be g_freed when done.
 */
static olong* medianPoln (ObitUVDesc *inDesc, ObitUVDesc *outDesc)
{
  olong *out=NULL;
  olong i, jincs, jincf, jincif, nif, pols[] = {1,2,4,8};
  olong ichan, iif, ipoln;

  out = g_malloc0(inDesc->ncorr*sizeof(olong));  /* Allocate output array */

  /* Increments */
  jincs = inDesc->incs / inDesc->inaxes[0];
  jincf = inDesc->incf / inDesc->inaxes[0];
  if (inDesc->jlocif>0) {
    jincif = inDesc->incif / inDesc->inaxes[0];
    nif    = inDesc->inaxes[inDesc->jlocif];
  }  else {
    jincif = 1;
    nif    = 1;
  }
  
  /* Loop over IF */
  for (iif=0; iif<nif; iif++) {
    /* Loop over channels */
    for (ichan=0; ichan<inDesc->inaxes[inDesc->jlocf]; ichan++) {
      /* Loop over poln */
      for (ipoln=0; ipoln<inDesc->inaxes[inDesc->jlocs]; ipoln++) {
	i = ichan*jincf + iif*jincif + ipoln*jincs;
	if (inDesc->crval[inDesc->jlocs]==1.0) out[i] = 15;  /* Stokes I */
	else out[i] = pols[ipoln];
      }
    }
  }					       
return out;
} /* end medianPoln */

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

/** 
 * Determine MedianDev function in a possibly threaded fashion
 * Arguments are given in the structure passed as arg
 * \param arg  Pointer to UVMednEditFuncArg argument with elements
 * \li thread  Thread with restart queue
 * \li ithread Thread number, >0 -> no threading 
 * \li first   First (1-rel) baseline to process this thread 
 * \li last    Highest (1-rel) baseline to process this thread 
 * \li amps    Circular amplitude storage (baseline,correlator,time)
 * \li itime   Time slice of interest 
 * \li times   Times of data samples, <-1000 => ignore
 * \li numBL   Number of baselines
 * \li numCorr Number of correlations
 * \li numTime Number of allocated time slices in amps
 * \li ntime   Number of actual time slices
 * \li alpha   Controls averaging
 *               0 -> 1 = pure boxcar -> pure MWF (alpha of the center data samples 
 *               are ignored and the rest averaged). 
 * \li devs    [out] Deviations in sigma per baseline/correlator, fblank => not determined 
 * \li work    Work array at least the size of ntime
 * \li number  [out] number of baselines/correlations with valid data 
 */
static gpointer ThreadMedianDev (gpointer arg)
{
  /* Get arguments from structure */
  UVMednEditFuncArg *largs = (UVMednEditFuncArg*)arg;
  olong loBL      = largs->first-1;
  olong hiBL      = largs->last;
  ofloat *amps    = largs->amps;
  olong itime     = largs->itime;  
  ofloat* times   = largs->times;   
  olong numCorr   = largs->numCorr; 
  olong numTime   = largs->numTime; 
  olong ntime     = largs->ntime;  
  ofloat alpha    = largs->alpha;   
  ofloat* devs    = largs->devs;   
  ofloat* work    = largs->work;

  ollong out = 0;
  olong  iBL, icorr, it, count, jndx, indx;
  ofloat delta, sigma, level, fblank = ObitMagicF();

  /* Loop over baselines */
  for (iBL=loBL; iBL<hiBL; iBL++) {
    /* Loop over correlator */
    for (icorr=0; icorr<numCorr; icorr++) {

      /* Copy valid data to work */
      count = 0;
      for (it=0; it<ntime; it++) {
	if (times[it]<-1000.0) continue; /* ignore blanked times */
	jndx = numTime * (icorr + iBL*numCorr) + it;
	if (amps[jndx]!=fblank) work[count++] = amps[jndx];
      } /* end loop over time */

      indx = icorr + iBL*numCorr;
      jndx = numTime * (icorr + iBL*numCorr) + itime;

      /* any data? */
      if ((count<=3) || (amps[jndx]==fblank)) {  /* no */
	devs[indx] = fblank;
      } else {         /* yes */
	level = MedianLevel (count, work, alpha);
	delta = fabs(amps[jndx] - level);
	/* This sigma will be an underestimate - double */
	sigma = 2*MedianSigma (count, work, level, alpha);
	/* Don't go overboard - min 0.1% */
	if (level!=fblank) sigma = MAX (0.001*level, sigma);
	if ((sigma>0.0) && (level!=fblank)) {
	  devs[indx] = delta/sigma;
	  out++;   /* Count valid data */
	} else devs[indx] = fblank;
      }
    } /* end loop over correlator  */
  } /* end loop over baseline */

  largs->number = out;  /* Save number */

  /* Indicate completion */
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadMedianDev */

