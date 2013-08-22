/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2013                                          */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include <gsl/gsl_multifit.h>
#include <stdlib.h>
#include "ObitOTFGetSoln.h"
#include "ObitTableOTFFlag.h"
#include "ObitOTFUtil.h"
#include "ObitTimeFilter.h"
#include "ObitPennArrayAtmFit.h"
#include "ObitPlot.h"
#include "ObitUtil.h"
#include "ObitTableOTFTargetUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFGetSoln.c
 * ObitOTFGetSoln utility module calibration function definitions.
 */

#define MAXSAMPSCAN 5000   /* Maximum number of samples in a scan */

/*---------------Private function prototypes----------------*/
/** Private: Solve for atmospheric offset calibration */
static void 
ObitOTFGetSolnSolve (ObitOTFArrayGeom *geom, ObitOTFDesc *desc, 
		     ObitPennArrayAtmFit *fitter, ofloat *rec, ObitTableOTFSolnRow *row);

/** Private: Solve for gain calibration */
static void 
ObitOTFGetSolnGainSolve (olong nsample, ofloat *time, ofloat *cal, ofloat calJy, ofloat *data, 
			 ofloat minrms, ofloat lastgood, olong detect, ObitTableOTFSolnRow *row);

/** Private: Solve for instrumental calibration */
static void 
ObitOTFGetInstSolve (ObitOTFArrayGeom *geom, ObitOTFDesc *desc, 
		     olong nDet, olong nTime, ofloat **data, ofloat *iCal, 
		     ObitTableOTFSolnRow *row);

/** Private: Get average around median value of an array */
static ofloat GetSolnMedianAvg (ofloat *array, olong incs, olong n);

/** Private: Fit polynomial */
static void  FitBLPoly (ofloat *poly, olong order, ofloat *x, 
			ofloat *y, ofloat *wt, olong n);

/** Private: Fit polynomial */
static void  FitMBBLPoly (ofloat solint, olong *npoly, ofloat **tpoly, ofloat **poly, 
			  ofloat *offset, olong ndetect, olong maxdata, 
			  ofloat *x, ofloat *y, ofloat *wt, olong n);
/** Private: Reject outlyers */
static void  FitMBBLOut (olong npoly, ofloat *tpoly, ofloat *poly, ofloat *offset, 
			 olong ndetect, olong maxdata, ofloat *x, ofloat *y, 
			 ofloat *wt, olong n, ofloat sigma);

/** Private: Plot multibeam baseline data and model */
static void PlotMBBL (olong npoly, ofloat *tpoly, ofloat *poly, ofloat *offset, 
		      olong ndetect, olong maxdata, 
		      ofloat *x, ofloat *y, ofloat *wt, olong n,
		      olong plotDetect, ofloat t0, ObitErr *err);

/** Private: Flatten curves, determine cal and Weights */
static void doPARCalWeight (olong nDet, olong nTime, ofloat *iCal, ofloat **accum, 
			    ofloat *time, ofloat *lastCal, ofloat *lastWeight, 
			    gboolean fitWate);

/** Private: qsort ofloat comparison */
static int compare_gfloat  (const void* arg1,  const void* arg2);

/** return opacity */
static ofloat getTau0 (ofloat time, gint32 taudim[], ofloat *taudata);

/** return RMS of an array (with blanking) */
static ofloat RMSValue (ofloat *array, olong n);

/** return ID of "Blank" Scan  */
static ofloat FindBlankID (ObitOTF *inOTF, ObitErr *err);

/*----------------------Public functions---------------------------*/
/**
 * Determine offset calibration for an OTF residual for multibeam
 * systems using an Atmospheric model across the array of detectors.
 * Calibration parameters are on the inOTF info member.
 * \li "solInt"    OBIT_float (1,1,1) Solution interval in days [def 1 sec].
 * \li "Tau0"      OBIT_float (?,?,1) Zenith opacity in nepers [def 0].
 *                 Can be passed as either a constant scalar or
 *                 an array of (time, tau0) pairs.
 * \li "minEl"     OBIT_float (1,1,1) Minimum elevation allowed (deg)
 * \li "Clip"      OBIT_float (1,1,1) data outside of range +/- Clip are replaced by
 *                 + or - Clip. [Def 1.0e20]
 * \li "plotTime"  OBIT_boolean (1,1,1) If True plot filtered time serieses
 * \li "plotFreq"  OBIT_boolean (1,1,1) If True plot filtered power spectra
 * \param inOTF    Input OTF data. 
 * \param outOTF   OTF with which the output  OTFSoln is to be associated
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created OTFSoln object which is 
 * associated with outOTF.
 */
ObitTableOTFSoln* ObitOTFGetSolnCal (ObitOTF *inOTF, ObitOTF *outOTF, ObitErr *err)
{
  ObitTableOTFSoln *outSoln=NULL;
  ObitTableOTFSolnRow *row=NULL;
  ObitPennArrayAtmFit *fitter;
  ObitTimeFilter *filter = NULL;
  gint32 taudim[MAXINFOELEMDIM] = {1,1,1,1,1}, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ObitIOAccess access;
  ofloat *rec=NULL, solInt, t0, fblank = ObitMagicF();
  ofloat *tau0, tMult, airmass;
  ofloat minEl, el, clip, lastScan=-1.0, lastTarget=-1.0;
  ofloat *time, dTime, *data, coef[20];
  ofloat totalTime, samp, parms[10];
  olong iRow, ver, i, j, k, ibuf, lrec, ncoef, nsample;
  olong  npoly, ndetect, incdatawt;
  gboolean flag, doCalSelect, someOK, someData, allBad;
  gboolean plotTime, plotFreq;
  ObitIOCode retCode;
  gchar *tname;
  gchar *routine = "ObitOTFGetSolnCal";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outSoln;
  g_assert (ObitOTFIsA(inOTF));

  /* Calibration wanted? */ 
  doCalSelect = FALSE;
  ObitInfoListGetTest(inOTF->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

  /* open OTF data if not already open */
  if ((inOTF->myStatus==OBIT_Inactive) || (inOTF->myStatus==OBIT_Defined)) {
    retCode = ObitOTFOpen (inOTF, access, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
  }

  /* Create output */
  tname = g_strconcat ("Calibration for: ",inOTF->name, NULL);
  ver = 0;
  ncoef = 1;  /* do first order atmospheric model (0=1, 1=3, 2nd = 6, 3rd=10) */
  npoly = ncoef;
  outSoln = newObitTableOTFSolnValue(tname, (ObitData*)outOTF, &ver, OBIT_IO_WriteOnly,  
				   inOTF->geom->numberDetect, npoly, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Create work arrays */
  ndetect = inOTF->geom->numberDetect;  /* number of detectors */
  lrec = inOTF->myDesc->lrec;
  time = g_malloc0(MAXSAMPSCAN*sizeof(ofloat));
  data = g_malloc0(ncoef*MAXSAMPSCAN*sizeof(ofloat));
  nsample  = 0;
  t0       = -1.0e20;
  lastScan = -1000.0;
  incdatawt = inOTF->myDesc->incdatawt; /* increment in data-wt axis */

  /* Get parameters for calibration */
  /* Solution interval default 10 sec */
  solInt = 10.0 / 86400.0;
  ObitInfoListGetTest(inOTF->info, "solInt", &type, dim, (gpointer*)&solInt);

  /* minimum allowed elevation for solution */
  minEl = 1.0;
  ObitInfoListGetTest(inOTF->info, "minEl",  &type, dim, (gpointer*)&minEl);

  /* Clip value */
  clip = 1.0e20;
  ObitInfoListGetTest(inOTF->info, "Clip",  &type, dim, (gpointer*)&clip);
  if (clip<1.0e19) 
    Obit_log_error(err, OBIT_InfoErr, "%s: Clipping residuals at %f", routine,clip);

  /* Opacity */
  tau0 = NULL;
  ObitInfoListGetP(inOTF->info, "Tau0",   &type, taudim, (gpointer*)&tau0);

  /* Plotting? */
  plotTime = FALSE;
  ObitInfoListGetTest(inOTF->info, "plotTime",  &type, dim, (gpointer*)&plotTime);
  plotFreq = FALSE;
  ObitInfoListGetTest(inOTF->info, "plotFreq",  &type, dim, (gpointer*)&plotFreq);

  /* Open output table */
  if ((ObitTableOTFSolnOpen (outSoln, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "%s ERROR opening input OTFSoln table", routine);
    return outSoln;
  }

  /* Create Row */
  row = newObitTableOTFSolnRow (outSoln);

  /* Attach row to output buffer */
  ObitTableOTFSolnSetRow (outSoln, row, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Initialize solution row */
  row->dAz = 0.0;  /* no pointing corrections */
  row->dEl = 0.0; /* no pointing corrections */
  row->Target = 0;
  for (j=0; j<ndetect; j++) row->mult[j] = 1.0;
  for (j=0; j<ndetect; j++) row->wt[j]   = 1.0;
  for (j=0; j<ndetect; j++) row->cal[j]  = 0.0;
  for (j=0; j<ndetect; j++) row->add[j]  = 0.0;
  for (i=0; i<npoly; i++)   row->poly[i] = 0.0;
  flag     = FALSE;
  someOK   = FALSE;
  someData = FALSE;

  /* Create fitter for atmospheric model */
  fitter = ObitPennArrayAtmFitValue ("Atm Fitter", inOTF->geom, ncoef, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
  
  /* loop over input data */
  retCode = OBIT_IO_OK;
  while (retCode == OBIT_IO_OK) {
    
    /* read buffer */
    if (doCalSelect) retCode = ObitOTFReadSelect (inOTF, NULL, err);
    else retCode = ObitOTFRead (inOTF, NULL, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
    if (retCode==OBIT_IO_EOF) break; /* done? */
    
    /* Record pointer */
    rec = inOTF->buffer;
    
    if (lastScan<0.0) {
      lastScan   = rec[inOTF->myDesc->ilocscan];  /* Which scan number */
      lastTarget = rec[inOTF->myDesc->iloctar];   /* Which target number */
    }
    
    /* Loop over buffer */
    for (ibuf=0; ibuf<inOTF->myDesc->numRecBuff; ibuf++) {
      
      /* Scan read? If so, compute solution and write. Also if work arrays full */
      if (( rec[inOTF->myDesc->ilocscan] != lastScan) || /* new scan */
	  (nsample>=MAXSAMPSCAN)) {   /* or exceed limit on number of samples */
	
	/* Any good data to filter? */
	if (someOK) {
	  
	  /* Filter data */
	  /* clip data */
	  if (clip<1.0e19) {
	    allBad = TRUE;
	    for (j=0; j<npoly; j++) {
	      for (k=0; k<nsample; k++) {
		if (data[j*MAXSAMPSCAN+k]!=fblank) {
		  if (data[j*MAXSAMPSCAN+k] > clip) data[j*MAXSAMPSCAN+k] = fblank;
		  if (data[j*MAXSAMPSCAN+k] <-clip) data[j*MAXSAMPSCAN+k] = fblank;
		  if (data[j*MAXSAMPSCAN+k]!=fblank) allBad = FALSE;
							}
	      }
	    }
	    /* Warn and flag if all clipped */
	    if (allBad) {
	      someOK = FALSE;
	      Obit_log_error(err, OBIT_InfoWarn, 
			     "Warning: All data in a scan clipped");
	    }
	  } /* End if clipping */
	  
	  /* Create filter as needed */
	  if (filter==NULL) filter = newObitTimeFilter("Cal Filter", nsample, npoly);
	  
	  /* Copy data to filter */
	  dTime = (time[nsample-1] - time[0]) / nsample;  /* Time increment */
	  for (j=0; j<npoly; j++) {
	    ObitTimeFilterGridTime (filter, j, dTime, nsample, time, &data[j*MAXSAMPSCAN]);
	  }

	  /* Plot? */
	  if (plotTime) ObitTimeFilterPlotTime (filter, 0, "Before filter", err);
	  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
	  
	  /* Transform to frequency */
	  ObitTimeFilter2Freq (filter);
	  
	  /* Apply filter */
	  totalTime = (time[nsample-1] - time[0]);  /* Total time in scan */
	  samp = totalTime / nsample;               /* sampling interval */
	  parms[0] = 2.0*samp / MAX(solInt,samp);   /* Cutoff frequency in units of highest */

	  parms[0] = 1.0/ (solInt*86400.0);         /* Cutoff frequency in Hz */
	  ObitTimeFilterDoFilter (filter, -1, OBIT_TimeFilter_LowPass, parms, err);
	  if (err->error) {
	    Obit_log_error(err, OBIT_InfoErr, "%s: Total time %f sample time %f No. %d solInt %f", 
			   routine,totalTime, samp, nsample, solInt);
	    Obit_traceback_val (err, routine, inOTF->name, outSoln);
	  }
	  
	  /* Plot? */
	  if (plotFreq) ObitTimeFilterPlotPower (filter, 0, "Power", err);
	  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
	  
	  /* FT back to time */
	  ObitTimeFilter2Time (filter);
	
	  /* Plot? */
	  if (plotTime) ObitTimeFilterPlotTime (filter, 0, "After filter", err);
	  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
	  
	  /* Copy data back to time series */
	  for (j=0; j<npoly; j++) {
	    ObitTimeFilterUngridTime (filter, j, nsample, time, &data[j*MAXSAMPSCAN]);
	  }
	} /* end of filter good data */
      
      /* Loop over solutions in scan */
      k = 0;
      t0 = time[0];
      while (k<nsample) {
	
	/* Use values in filtered time series separated by solInt/4 
	   to over Nyquist sample highest frequency */
	/* Find next time */
	while ((time[k]<t0) && (k<nsample)) {k++;}
	k = MIN (k, (nsample-1));
	
	/* Set descriptive info on Row */
	row->Time  = time[k];      /* time */
	row->TimeI = 0.25*solInt;  /* time interval*/
	row->Target = (oint)(lastTarget+0.5);

	/* Opacity correction */
	el = ObitOTFArrayGeomElev (inOTF->geom, row->Time, 
				 rec[inOTF->myDesc->ilocra], 
				 rec[inOTF->myDesc->ilocdec]);
	if (el > 0.1) 
	  airmass = 1.0 / cos (DG2RAD * (90.0 - el));
	else
	  airmass = 10.0;
	tMult = exp (getTau0(row->Time,taudim,tau0) * airmass); /* multiplicative term */
	for (j=0; j<ndetect; j++) row->mult[j] = tMult;
	
	/* Write Soln table */
	iRow = -1;
	/* Set poly filtered values */
	if (someOK) {  /* Some data OK */
	  for (j=0; j<npoly; j++)  row->poly[j] = data[j*MAXSAMPSCAN+k];
	} else { /* no good data - blank */
	  for (j=0; j<npoly; j++)  row->poly[j] = fblank;
	}
	if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	     != OBIT_IO_OK) || (err->error>0)) { 
	  Obit_log_error(err, OBIT_Error, "%s ERROR writing OTFSoln Table file", routine);
	  return outSoln;
	}
	
	t0 += 0.25*solInt; /* time of next entry in Soln table */
	
	/* done? */
	if (k>=(nsample-1)) break;
	
      } /* end loop over times in scan */
      
	/* initialize accumulators */
      t0 = rec[inOTF->myDesc->iloct];
      nsample = 0;
      lastScan = -1000.0;
      someOK = FALSE;
      
      } /* end do solution */
      
      /* Save sata */
      /* Is data below the elevation limit? */
      el = ObitOTFArrayGeomElev (inOTF->geom, rec[inOTF->myDesc->iloct], 
				 rec[inOTF->myDesc->ilocra], 
				 rec[inOTF->myDesc->ilocdec]);
      flag = (el < minEl);

      /* Fit model per data sample */
      ObitPennArrayAtmFitFit (fitter, &rec[inOTF->myDesc->ilocdata], incdatawt, coef);
     
      /* accumulate */
      time[nsample] = rec[inOTF->myDesc->iloct];
      if (!flag) { /* OK */
	for (j=0; j<npoly; j++ ) {
	  data[nsample+j*MAXSAMPSCAN] = coef[j];
	  someOK = someOK || (data[nsample+j*MAXSAMPSCAN]!=fblank);
	}
      } else { /* too low el - blank */
	for (j=0; j<npoly; j++ ) data[nsample+j*MAXSAMPSCAN] = fblank;
      }
      nsample++; /* how many samples in buffer */
      someData = someData || someOK;  /* Any valid data? */
      
      if (lastScan<0.0) {
	lastScan   = rec[inOTF->myDesc->ilocscan];  /* Which scan number */
	lastTarget = rec[inOTF->myDesc->iloctar];   /* Which target number */
      }
      rec += inOTF->myDesc->lrec; /* Update data record pointer */
      
    } /* end loop over buffer load */
  } /* end loop reading data */
  
  /* Do last scan still in arrays*/
  if (nsample>1) {

    /* Any good data to filter? */
    if (someOK) {

	  /* Filter data */
	  /* clip data */
	  if (clip<1.0e19) {
	    allBad = TRUE;
	    for (j=0; j<npoly; j++) {
	      for (k=0; k<nsample; k++) {
		if (data[j*MAXSAMPSCAN+k]!=fblank) {
		  if (data[j*MAXSAMPSCAN+k] > clip) data[j*MAXSAMPSCAN+k] = fblank;
		  if (data[j*MAXSAMPSCAN+k] <-clip) data[j*MAXSAMPSCAN+k] = fblank;
		  if (data[j*MAXSAMPSCAN+k]!=fblank) allBad = FALSE;
							}
	      }
	    }
	    /* Warn and flag if all clipped */
	    if (allBad) {
	      someOK = FALSE;
	      Obit_log_error(err, OBIT_InfoWarn, 
			     "Warning: All data in a scan clipped");
	    }
	  } /* End if clipping */
	  
	  /* Create filter as needed */
	  if (filter==NULL) filter = newObitTimeFilter("Cal Filter", nsample, npoly);
	  
	  /* Copy data to filter */
	  dTime = (time[nsample-1] - time[0]) / nsample;  /* Time increment */
	  for (j=0; j<npoly; j++) {
	    ObitTimeFilterGridTime (filter, j, dTime, nsample, time, &data[j*MAXSAMPSCAN]);
	  }
	  
	  /* Transform to frequency */
	  ObitTimeFilter2Freq (filter);
	  
	  /* Apply filter */
	  totalTime = (time[nsample-1] - time[0]);  /* Total time in scan */
	  samp = totalTime / nsample;               /* sampling interval */
	  parms[0] = 2.0*samp / MAX(solInt,samp);   /* Cutoff frequency in units of highest */

	  parms[0] = 1.0/ (solInt*86400.0);         /* Cutoff frequency in Hz */
	  ObitTimeFilterDoFilter (filter, -1, OBIT_TimeFilter_LowPass, parms, err);
	  if (err->error) {
	    Obit_log_error(err, OBIT_InfoErr, "%s: Total time %f sample time %f No. %d solInt %f", 
			   routine,totalTime, samp, nsample, solInt);
	    Obit_traceback_val (err, routine, inOTF->name, outSoln);
	  }
	  
	  /* Plot? */
	  if (plotTime) ObitTimeFilterPlotTime (filter, 0, "Before filter", err);
	  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
	  
	  /* FT back to time */
	  ObitTimeFilter2Time (filter);
	
	  /* Time Plot? */
	  if (plotTime) ObitTimeFilterPlotTime (filter, 0, "After filter", err);
	  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
	  
	  /* Power spectrum plot? */
	  if (plotTime) ObitTimeFilterPlotTime (filter, 0, "Power", err);
	  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
	  
	  /* Copy data back to time series */
	  for (j=0; j<npoly; j++) {
	    ObitTimeFilterUngridTime (filter, j, nsample, time, &data[j*MAXSAMPSCAN]);
	  }
	} /* end of filter good data */
      
      /* Loop over solutions in scan */
      k = 0;
      t0 = time[0];
      while (k<nsample) {
	
	/* Use values in filtered time series separated by solInt/4 
	   to over Nyquist sample highest frequency */
	/* Find next time */
	while ((time[k]<t0) && (k<nsample)) {k++;}
	k = MIN (k, (nsample-1));
	
	/* Set descriptive info on Row */
	row->Time  = time[k];      /* time */
	row->TimeI = 0.25*solInt;  /* time interval*/
	row->Target = (oint)(lastTarget+0.5);
	
	/* Opacity correction */
	el = ObitOTFArrayGeomElev (inOTF->geom, row->Time, 
				 rec[inOTF->myDesc->ilocra], 
				 rec[inOTF->myDesc->ilocdec]);
	if (el > 0.1) 
	  airmass = 1.0 / cos (DG2RAD * (90.0 - el));
	else
	  airmass = 10.0;
	tMult = exp (getTau0(row->Time,taudim,tau0) * airmass); /* multiplicative term */
	for (j=0; j<ndetect; j++) row->mult[j] = tMult;
	
	/* Write Soln table */
	iRow = -1;
	/* Set poly filtered values */
	if (someOK) {  /* Some data OK */
	  for (j=0; j<npoly; j++)  row->poly[j] = data[j*MAXSAMPSCAN+k];
	} else { /* no good data - blank */
	  for (j=0; j<npoly; j++)  row->poly[j] = fblank;
	}
	if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	     != OBIT_IO_OK) || (err->error>0)) { 
	  Obit_log_error(err, OBIT_Error, "%s ERROR writing OTFSoln Table file", routine);
	  return outSoln;
	}
	
	t0 += 0.25*solInt; /* time of next entry in Soln table */
	
	/* done? */
	if (k>=(nsample-1)) break;
	
      } /* end loop over times in scan */
      
  } /* end finish up data in arrays */

  
  /* Close output cal table */
  if ((ObitTableOTFSolnClose (outSoln, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s ERROR closing output OTFSoln Table file", routine);
    return outSoln;
  }
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Give warning if no data selected */
  if (!someData) Obit_log_error(err, OBIT_InfoWarn, 
				"%s: Warning: NO data selected", routine);

  /* Cleanup */
  fitter = ObitPennArrayAtmFitUnref (fitter);
  filter = ObitTimeFilterUnref(filter);
  row = ObitTableOTFSolnUnref(row);
  if (time) g_free(time);
  if (data) g_free(data);

  return outSoln;
} /* end ObitOTFGetSolnCal */

/**
 * Determine gain calibration for an OTF from a residual data set.
 * Gain Calibration is based on the "Cal" values in the data.
 * Average value of the noise Cal is computed and entered in the data. (Gain only)
 * Additive values are determined from the median values of the residual data.
 * Solution type controlled by calType 
 * Calibration parameters are on the inOTF info member.
 * \li "solInt"    OBIT_float (1,1,1) Solution interval in days [def 10 sec].
 *                 This should not exceed 1000 samples.  Solutions will be truncated
 *                 at this limit.
 * \li "minRMS"    OBIT_float (1,1,1) minimum allowable fractional solution RMS [def 0.1].
 *                 bad solutions are replaced with pervious good value. [Gain soln]
 * \li "calJy"     OBIT_float (*,1,1) Calibrator value in Jy per detector [Gain soln] [def 1.0] .
 * \li "minEl"     OBIT_float (1,1,1) Minimum elevation allowed (deg)
 * \li "calType"   OBIT_string (*,1,1) Calibration type desired
 *                 "Gain" => Gain (multiplicative, cal) solution only
 *                 "Offset" => Offset (additive) solution only
 *                 "GainOffset" both gain and offset calibration (probably bad idea).
 *                 anything else or absent => Gain only.
 * \param inOTF    Input OTF data. 
 * \param outOTF   OTF with which the output  OTFSoln is to be associated
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created OTFSoln object which is 
 * associated with inOTF.
 */
ObitTableOTFSoln* ObitOTFGetSolnGain (ObitOTF *inOTF, ObitOTF *outOTF, ObitErr *err)
{
#define MAXSAMPLE 1000   /* Maximum number of samples in an integration */
  ObitTableOTFSoln *outSoln=NULL;
  ObitTableOTFSolnRow *row=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ObitOTFDesc *desc=inOTF->myDesc;
  ofloat lastScan=-1.0, lastTarget=-1.0;
  ofloat *rec=NULL, solInt, minEl, el, t0, tEnd, value, fblank = ObitMagicF();
  ofloat *time, *cal, *data, *lastgood, sumTime, minrms, *calJy=NULL;
  ofloat *value1=NULL, *value2=NULL;
  olong nfeed, nstok, nchan, nscal, ifeed, istok, incs, incfeed, ilocdata;
  olong iRow, ver, i, j, lrec, nsample;
  olong  npoly, ndetect, incdatawt;
  ObitIOCode retCode;
  gboolean flag, someData, doGain, doOffset;
  gchar *tname;
  gchar calType[50];
  gchar *routine = "ObitOTFGetSolnGain";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outSoln;
  g_assert (ObitOTFIsA(inOTF));

  /* open OTF data to fully instantiate if not already open */
  if ((inOTF->myStatus==OBIT_Inactive) || (inOTF->myStatus==OBIT_Defined)) {
    retCode = ObitOTFOpen (inOTF, OBIT_IO_ReadWrite, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
  }

  /* Create output */
  tname = g_strconcat ("Calibration for: ",inOTF->name, NULL);
  ver = 0;
  npoly = 1;
  outSoln = newObitTableOTFSolnValue(tname, (ObitData*)outOTF, &ver, OBIT_IO_WriteOnly,  
				   inOTF->geom->numberDetect, npoly, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* How many detectors? */
  if (desc->OTFType==OBIT_GBTOTF_VEGAS) {
    incs      = desc->incs;
    incfeed   = desc->incfeed;
    ilocdata  = desc->ilocdata;
    nfeed     = desc->inaxes[desc->jlocfeed];
    nstok     = desc->inaxes[desc->jlocs];
    nscal     = MIN (2, nstok);              /* Number of stokes in cal table */
    nchan     = desc->inaxes[desc->jlocf];
    ndetect   = nscal*nfeed;
    /* Better be frequency averaged */
    Obit_retval_if_fail ((nchan==1), err, outSoln,
			 "%s VEGAS data must be frequency averaged",  
			 routine);  
  } else { /* Non VEGAS */
    ndetect = inOTF->geom->numberDetect;  /* number of detectors */
  }
  /* Create work arrays */
  lrec     = inOTF->myDesc->lrec;
  lastgood = g_malloc0(ndetect*sizeof(ofloat));
  time     = g_malloc0(MAXSAMPLE*sizeof(ofloat));
  cal      = g_malloc0(MAXSAMPLE*sizeof(ofloat));
  data     = g_malloc0(ndetect*MAXSAMPLE*sizeof(ofloat));
  calJy    = g_malloc0(ndetect*sizeof(ofloat));
  value1   = g_malloc0(ndetect*sizeof(ofloat));
  value2   = g_malloc0(ndetect*sizeof(ofloat));
  nsample  = 0;
  t0       = -1.0e20;
  tEnd     = -1.0e20;
  sumTime  = 0.0;
  lastScan = -1000;
  incdatawt = inOTF->myDesc->incdatawt; /* increment in data-wt axis */
  someData = FALSE;

  /* Get parameters for calibration */
  /* Solution interval default 10 sec */
  solInt = 10.0 / 86400.0;
  ObitInfoListGetTest(inOTF->info, "solInt", &type, dim, (gpointer*)&solInt);

  /* minimum RMS in solution */
  minrms = 0.1;
  ObitInfoListGetTest(inOTF->info, "minRMS", &type, dim, (gpointer*)&minrms);

  /* minimum allowed elevation for solution */
  minEl = 1.0;
  ObitInfoListGetTest(inOTF->info, "minEl",  &type, dim, (gpointer*)&minEl);

  /* cal value in Jy */
  for (j=0; j<ndetect; j++) calJy[j] = 1.0;
  ObitInfoListGetTest(inOTF->info, "calJy",  &type, dim, (gpointer*)calJy);
 
  /* calibration type */
  for (i=0; i<50; i++) calType[i] = 0;
  strcpy (calType, "Gain");
  ObitInfoListGetTest(inOTF->info, "calType", &type, dim, calType);
  /* What calibrations are desired? */
  doGain = FALSE; doOffset = FALSE;
  if (!strncmp(calType, "Gain",4)) doGain = TRUE;
  else if (!strncmp(calType, "Offset",6)) doOffset = TRUE;
  else if (!strncmp(calType, "GainOffset",10)) {doGain=TRUE; doOffset = TRUE;}
  else doGain = TRUE;

  /* Open table */
  if ((ObitTableOTFSolnOpen (outSoln, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input OTFSoln table");
    return outSoln;
  }

  /* Create Row */
  row = newObitTableOTFSolnRow (outSoln);

  /* Attach row to output buffer */
  ObitTableOTFSolnSetRow (outSoln, row, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Initialize solution row */
  row->dAz = 0.0;  /* no pointing corrections */
  row->dEl = 0.0; /* no pointing corrections */
  row->Target = 0;
  for (j=0; j<ndetect; j++) row->mult[j] = 1.0;
  for (j=0; j<ndetect; j++) row->wt[j]   = 1.0;
  for (j=0; j<ndetect; j++) row->cal[j]  = 0.0;
  for (j=0; j<ndetect; j++) row->add[j]  = 0.0;
  for (j=0; j<ndetect; j++) lastgood[j]  = fblank;
  for (j=0; j<npoly; j++)   row->poly[j] = 0.0;
  flag = FALSE;
  rec = inOTF->buffer;

  /* loop calibrating data */
  retCode = OBIT_IO_OK;
  while (retCode == OBIT_IO_OK) {

    /* read buffer */
    retCode = ObitOTFRead (inOTF, NULL, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
    if (retCode==OBIT_IO_EOF) break; /* done? */

    /* Record pointer */
    rec = inOTF->buffer;
  
    /* First time */
    if (t0<-1.0e10) {
      t0 = rec[inOTF->myDesc->iloct]; 
      tEnd = t0+1.0E-20;
      lastScan   = rec[inOTF->myDesc->ilocscan]; /* Which scan number */
      lastTarget = rec[inOTF->myDesc->iloctar];  /* Which target number */
    }

    /* Loop over buffer */
    for (i=0; i<inOTF->myDesc->numRecBuff; i++) {

      /* Solution interval finished? If so, compute solution and write.
	 also if work arrays full */
      if ((rec[inOTF->myDesc->iloct] > (t0+solInt)) || 
	  /* or end of scan */
	  ( rec[inOTF->myDesc->ilocscan] != lastScan) ||
	  /* or exceed limit on number of samples */
	  (nsample>=MAXSAMPLE)) {
	
	/* Set descriptive info on Row */
	row->Target = (oint)(lastTarget+0.5);
	if (nsample>0) {
	  row->Time  = sumTime/nsample;  /* time */
	  row->TimeI = 2.0 * (row->Time - t0);
	} else {
	  row->Time  = rec[inOTF->myDesc->iloct];
	  row->TimeI = 0.0;
	}
	
	/* Get calibration, per detector */
	for (j=0; j<ndetect; j++) {
	  row->add[j]  = 0.0;
	  row->mult[j] = 1.0;
	  row->wt[j]   = 1.0;
	  row->cal[j]  = 0.0;
	  if (doGain) {  /* gain solution */
	    ObitOTFGetSolnGainSolve (nsample, time, cal, calJy[j], &data[j*MAXSAMPLE], 
				     minrms, lastgood[j], j, row);
	    lastgood[j] = row->mult[j];  /* save last good value */
	  }
	  if (doOffset ) {  /* Offset solution */
	    if (nsample>1) {
	      /* Do it all at once */
	      value = GetSolnMedianAvg (&data[j*MAXSAMPLE], 1, nsample);
	      if (value!=fblank) value = -value;
	      value1[j] = value; /* for beginning */
	      value2[j] = value; /* for end */
	    } else { /* no data */
	      value1[j] = fblank;
	      value2[j] = fblank;
	    }
	  }
	}
	
	/* Write Soln table, write at beginning and end of solution */
	iRow = -1;
	/* Beginning */
	row->Time  = t0;
	/* Set offset solution values */
	if (doOffset) for (j=0; j<ndetect; j++)  row->add[j] = value1[j];
	if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	     != OBIT_IO_OK) || (err->error>0)) { 
	  Obit_log_error(err, OBIT_Error, "%s ERROR writing OTFSoln Table file", routine);
	  return outSoln;
	}
	/* End */
	row->Time  = tEnd;
	/* Set offset solution values */
	if (doOffset) for (j=0; j<ndetect; j++)  row->add[j] = value2[j];
	if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	     != OBIT_IO_OK) || (err->error>0)) { 
	  Obit_log_error(err, OBIT_Error, "%s ERROR writing OTFSoln Table file", routine);
	  return outSoln;
	}
	/* initialize accumulators */
	t0 = rec[inOTF->myDesc->iloct];
	tEnd = t0+1.0E-20;
	sumTime = 0.0;
	nsample = 0;
	
      } /* end do solution */
      
      /* Is data below the elevation limit? */
      el = ObitOTFArrayGeomElev (inOTF->geom, rec[inOTF->myDesc->iloct], 
				 rec[inOTF->myDesc->ilocra], 
				 rec[inOTF->myDesc->ilocdec]);
      flag = flag || (el < minEl);
      
      /* accumulate */
      if (!flag) {
	sumTime += rec[inOTF->myDesc->iloct];
	time[nsample] = rec[inOTF->myDesc->iloct];
	cal[nsample]  = rec[inOTF->myDesc->iloccal];
	/* VEGAS different use, parallel poln and feeds */
	j = 0;
	if (desc->OTFType==OBIT_GBTOTF_VEGAS) {
	  /* Loop over feed */
	  for (ifeed=0; ifeed<nfeed; ifeed++) {
	    /* Loop over parallel Stokes */
	    for (istok=0; istok<nscal; istok++) {
	       data[nsample+j*MAXSAMPLE] = 
		 rec[ilocdata + ifeed*incfeed + istok*incs];
	       j++;
	    } /* end Stokes loop */
	  } /* end feed loop */
	} else { /* Non VEGAS */
	  for (j=0; j<ndetect; j++ ) 
	    data[nsample+j*MAXSAMPLE] = rec[ilocdata+j*incdatawt];
	} /* end non VEGAS */
	nsample++;
	someData = TRUE;  /* Any valid data? */
      } /* end if not flagged */
      /* last time in accumulation */
      tEnd       = rec[inOTF->myDesc->iloct]+1.0E-20;
      lastScan   = rec[inOTF->myDesc->ilocscan]; /* Which scan number */
      lastTarget = rec[inOTF->myDesc->iloctar]; /* Which target number */
      rec += inOTF->myDesc->lrec; /* Data record pointer */
      
    } /* end loop over buffer load */
  } /* end loop reading data */
  
  /* Finish up any data in arrays */
  /* Set descriptive info on Row */
  row->Target = (oint)(lastTarget+0.5);
  if (nsample>0) {
    row->Time  = sumTime/nsample;  /* time */
    row->TimeI = 2.0 * (row->Time - t0);
  } else {
    row->Time  = rec[inOTF->myDesc->iloct];
    row->TimeI = 0.0;
  }
 
  /* Get calibration, per detector */
  for (j=0; j<ndetect; j++) {
    row->add[j]  = 0.0;
    row->mult[j] = 1.0;
    row->wt[j]   = 1.0;
    row->cal[j]  = 0.0;
    if (doGain) {  /* gain solution */
      ObitOTFGetSolnGainSolve (nsample, time, cal, calJy[j], &data[j*MAXSAMPLE], 
			       minrms, lastgood[j], j, row);
      lastgood[j] = row->mult[j];  /* save last good value */
    }
    if (doOffset ) {  /* Offset solution */
      if (nsample>1) {
	/* all at once */
	value = GetSolnMedianAvg (&data[j*MAXSAMPLE], 1, nsample);
	if (value!=fblank) value = -value;
	value1[j] = value;
	/* NO Do second half 
	   value = GetSolnMedianAvg (&data[j*MAXSAMPLE+nsample/2], incdatawt, nsample/2);
	   if (value!=fblank) value = -value;*/
	value2[j] = value;
      } else { /* no data */
	value1[j] = fblank;
	value2[j] = fblank;
      }
    }
  }

  /* Write Soln table, write at beginning and end of solution */
  iRow = -1;
  /* Beginning */
  row->Target = (oint)(lastTarget+0.5);
  row->Time  = t0;
  /* Set offset solution values */
  if (doOffset) for (j=0; j<ndetect; j++)  row->add[j]  = value1[j];
  if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "%s ERROR writing OTFSoln Table file", routine);
    return outSoln;
  }
  /* End */
  row->Time  = tEnd;
  /* Set offset solution values */
  if (doOffset) for (j=0; j<ndetect; j++)  row->add[j]  = value2[j];
  if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "%s ERROR writing OTFSoln Table file", routine);
    return outSoln;
  }
  /*  End final solution */

  
  /* Close cal table */
  if ((ObitTableOTFSolnClose (outSoln, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s ERROR closing output OTFSoln Table file", routine);
    return outSoln;
  }
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Give warning if no data selected */
  if (!someData) Obit_log_error(err, OBIT_InfoWarn, 
				"%s: Warning: NO data selected", routine);

  /* Cleanup */
  row = ObitTableOTFSolnUnref(row);
  if (time) g_free(time);
  if (cal) g_free(cal);
  if (data) g_free(data);

  return outSoln;
} /* end ObitOTFGetSolnGain */

/**
 * Determine offset calibration for an OTF by time filtering a residual data set.
 * The time series of each detector is filtered to remove structure on time 
 * scales longer than solInt. 
 * Scans in excess of 5000 samples will be broken into several.
 * Calibration parameters are on the inOTF info member.
 * \li "solInt"    OBIT_float (1,1,1) Solution interval in days [def 10 sec].
 *                 There will be 4 Soln table entries per solInt
 *                 This should not exceed 1000 samples.  Solutions will be truncated
 *                 at this limit.
 * \li "minEl"     OBIT_float (1,1,1) Minimum elevation allowed (deg)
 * \li "Clip"      OBIT_float (1,1,1) data outside of range +/- Clip are replaced by
 *                 + or - Clip. [Def 1.0e20]
 * \param inOTF    Input OTF data. Prior calibration applied if requested.
 * \param outOTF   OTF with which the output OTFSoln is to be associated
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created OTFSoln object which is 
 * associated with inOTF.
 */
ObitTableOTFSoln* ObitOTFGetSolnFilter (ObitOTF *inOTF, ObitOTF *outOTF, ObitErr *err)
{
  ObitTableOTFSoln *outSoln=NULL;
  ObitTableOTFSolnRow *row=NULL;
  ObitTimeFilter *filter = NULL;
  ObitOTFDesc *desc=inOTF->myDesc;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ObitIOAccess access;
  ofloat lastScan=-1.0, lastTarget=-1.0;
  ofloat *rec, solInt, minEl, el, t0, fblank = ObitMagicF();
  ofloat *time, *data, clip;
  ofloat totalTime, samp, parms[10];
  olong nfeed, nstok, nchan, nscal, ifeed, istok, incs, incfeed, ilocdata;
  olong iRow, ver, i, j, k, ibuf, lrec, nsample;
  olong  npoly, ndetect, incdatawt, nTime=0;
  ObitIOCode retCode;
  gboolean flag, someOK, someData, doCalSelect, done;
  gchar *tname;
  gchar *routine = "ObitOTFGetSolnFilter";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outSoln;
  g_assert (ObitOTFIsA(inOTF));

  /* Calibration wanted? */ 
  doCalSelect = FALSE;
  ObitInfoListGetTest(inOTF->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

  /* open OTF data if not already open */
  if ((inOTF->myStatus==OBIT_Inactive) || (inOTF->myStatus==OBIT_Defined)) {
    retCode = ObitOTFOpen (inOTF, access, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
  }

  /* Create output */
  tname = g_strconcat ("Calibration for: ",inOTF->name, NULL);
  ver = 0;
  npoly = 1;
  outSoln = newObitTableOTFSolnValue(tname, (ObitData*)outOTF, &ver, OBIT_IO_WriteOnly,  
				   inOTF->geom->numberDetect, npoly, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* How many detectors? */
  if (desc->OTFType==OBIT_GBTOTF_VEGAS) {
    incs      = desc->incs;
    incfeed   = desc->incfeed;
    ilocdata  = desc->ilocdata;
    nfeed     = desc->inaxes[desc->jlocfeed];
    nstok     = desc->inaxes[desc->jlocs];
    nscal     = MIN (2, nstok);              /* Number of stokes in cal table */
    nchan     = desc->inaxes[desc->jlocf];
    ndetect   = nscal*nfeed;
    /* Better be frequency averaged */
    Obit_retval_if_fail ((nchan==1), err, outSoln,
			 "%s VEGAS data must be frequency averaged",  
			 routine);  
  } else { /* Non VEGAS */
    ndetect = inOTF->geom->numberDetect;  /* number of detectors */
  }
  /* Create work arrays */
  lrec = inOTF->myDesc->lrec;
  time = g_malloc0(MAXSAMPSCAN*sizeof(ofloat));
  data = g_malloc0(ndetect*MAXSAMPSCAN*sizeof(ofloat));
  nsample = 0;
  t0   = -1.0e20;
  lastScan = -1000.0;
  incdatawt = inOTF->myDesc->incdatawt; /* increment in data-wt axis */

  /* Get parameters for calibration */
  /* Solution interval default 10 sec */
  solInt = 10.0 / 86400.0;
  ObitInfoListGetTest(inOTF->info, "solInt", &type, dim, (gpointer*)&solInt);

  /* minimum allowed elevation for solution */
  minEl = 1.0;
  ObitInfoListGetTest(inOTF->info, "minEl",  &type, dim, (gpointer*)&minEl);

  /* Clip value */
  clip = 1.0e20;
  ObitInfoListGetTest(inOTF->info, "Clip",  &type, dim, (gpointer*)&clip);
  if (clip<1.0e19) 
    Obit_log_error(err, OBIT_InfoErr, "%s: Clipping residuals at %f", routine,clip);

  /* Open output table */
  if ((ObitTableOTFSolnOpen (outSoln, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "%s ERROR opening input OTFSoln table", routine);
    return outSoln;
  }

  /* Create Row */
  row = newObitTableOTFSolnRow (outSoln);

  /* Attach row to output buffer */
  ObitTableOTFSolnSetRow (outSoln, row, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Initialize solution row */
  row->dAz = 0.0;  /* no pointing corrections */
  row->dEl = 0.0; /* no pointing corrections */
  row->Target = 0;
  for (j=0; j<ndetect; j++) row->mult[j] = 1.0;
  for (j=0; j<ndetect; j++) row->wt[j]   = 1.0;
  for (j=0; j<ndetect; j++) row->cal[j]  = 0.0;
  for (j=0; j<ndetect; j++) row->add[j]  = 0.0;
  for (i=0; i<npoly; i++)   row->poly[i] = 0.0;
  flag     = FALSE;
  someOK   = FALSE;
  someData = FALSE;

  /* loop over input data */
  retCode = OBIT_IO_OK;
  while (retCode == OBIT_IO_OK) {

    /* read buffer */
    if (doCalSelect) retCode = ObitOTFReadSelect (inOTF, NULL, err);
    else retCode = ObitOTFRead (inOTF, NULL, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
    done = retCode==OBIT_IO_EOF;
    if (done) break; /* done? */

    /* Record pointer */
    rec = inOTF->buffer;
  
    if (lastScan<0.0) {
      lastScan   = rec[inOTF->myDesc->ilocscan]; /* Which scan number */
      lastTarget = rec[inOTF->myDesc->iloctar];  /* Which target number */
    }

    /* Loop over buffer */
    for (ibuf=0; ibuf<inOTF->myDesc->numRecBuff; ibuf++) {

      /* Scan read? If so, compute solution and write. Also if work arrays full */
      if (( rec[inOTF->myDesc->ilocscan] != lastScan) || /* new scan */
	  (nsample>=MAXSAMPSCAN) || done) {   /* or exceed limit on number of samples, or done */
	
	/* Any good data to filter? */
	if (someOK) {
	  
	  /* Filter data */
	  /* Get optimum filter length - zero pad */
	  nTime = ObitFFTSuggestSize (nsample*2);
	  
	  /* Create or resize filter as needed */
	  if (filter==NULL) {
	    filter = newObitTimeFilter("Cal Filter", nTime, ndetect);
	  } else {
	    ObitTimeFilterResize (filter, nTime);
	  }
	  
	  /* copy data */
	  for (j=0; j<ndetect; j++) {
	    for (k=0; k<nsample; k++) filter->timeData[j][k] = data[j*MAXSAMPSCAN+k];
	  }
	  
	  /* clip data */
	  if (clip<1.0e19) {
	    for (j=0; j<ndetect; j++) {
	      for (k=0; k<nsample; k++) {
		if (filter->timeData[j][k]!=fblank) {
		  if (filter->timeData[j][k] > clip) filter->timeData[j][k] = clip;
		  if (filter->timeData[j][k] <-clip) filter->timeData[j][k] =-clip;
		}
	      }
	    }
	  }
	  
	  /* blank fill remainder of time series */
	  for (j=0; j<ndetect; j++) {
	    for (k=nsample; k<nTime; k++) filter->timeData[j][k] = fblank;
	  }
	  
	  /* Transform to frequency */
	  ObitTimeFilter2Freq (filter);
	  
	  /* Apply filter */
	  totalTime = (time[nsample-1] - time[0]);  /* Total time in scan */
	  samp = totalTime / nsample;               /* sampling interval */
	  parms[0] = 2.0 * samp / MAX(solInt,samp); /* Cutoff frequency in units of highest */
	  ObitTimeFilterFilter (filter, -1, OBIT_TimeFilter_LowPass, parms, err);
	  if (err->error) {
	    Obit_log_error(err, OBIT_InfoErr, "%s: Total time %f sample time %f No. %d solInt %f", 
			   routine,totalTime, samp, nsample, solInt);
	    Obit_traceback_val (err, routine, inOTF->name, outSoln);
	  }
	  
	  /* FT back to time */
	  ObitTimeFilter2Time (filter);
	  
	} /* end of filter good data */
	
	/* Loop over solutions in scan */
	k = 0;
	t0 = time[0];
	while (k<nsample) {
	  
	  /* Use values in filtered time series separated by solInt/4 
	     to over Nyquist sample highest frequency */
	  /* Find next time */
	  while ((time[k]<t0) && (k<nsample)) {k++;}
	  k = MIN (k, (nsample-1));
	  
	  /* Set descriptive info on Row */
	  row->Time  = time[k];      /* time */
	  row->TimeI = 0.25*solInt;  /* time interval*/
	  row->Target = (oint)(lastTarget+0.5);
	  
	  /* Write Soln table */
	  iRow = -1;
	  /* Set offset -filtered values */
	  if (someOK) {  /* Some data OK */
	    for (j=0; j<ndetect; j++)  row->add[j] = -filter->timeData[j][k];
	  } else { /* no good data - blank */
	    for (j=0; j<ndetect; j++)  row->add[j] = fblank;
	  }
	  if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	       != OBIT_IO_OK) || (err->error>0)) { 
	    Obit_log_error(err, OBIT_Error, "%s ERROR writing OTFSoln Table file", routine);
	    return outSoln;
	  }
	  
	  t0 += 0.25*solInt; /* time of next entry in Soln table */
	  
	  /* done? */
	  if (k>=(nsample-1)) break;
	  
	} /* end loop over times in scan */
	
	/* initialize accumulators */
	t0 = rec[inOTF->myDesc->iloct];
	nsample = 0;
	lastScan = -1000.0;
	someOK = FALSE;
	
      } /* end do solution */
      
      /* Is data below the elevation limit? */
      el = ObitOTFArrayGeomElev (inOTF->geom, rec[inOTF->myDesc->iloct], 
				 rec[inOTF->myDesc->ilocra], 
				 rec[inOTF->myDesc->ilocdec]);
      flag = (el < minEl);
      
      /* accumulate */
      time[nsample] = rec[inOTF->myDesc->iloct];
      if (!flag) { /* OK */
	/* VEGAS different use, parallel poln and feeds */
	j = 0;
	if (desc->OTFType==OBIT_GBTOTF_VEGAS) {
	  /* Loop over feed */
	  for (ifeed=0; ifeed<nfeed; ifeed++) {
	    /* Loop over parallel Stokes */
	    for (istok=0; istok<nscal; istok++) {
	      data[nsample+j*MAXSAMPSCAN] = 
		rec[ilocdata + ifeed*incfeed + istok*incs];
	      someOK = someOK || (data[nsample+j*MAXSAMPSCAN]!=fblank);
	      j++;
	    } /* end Stokes loop */
	  } /* end feed loop */
	} else { /* Non VEGAS */
	  for (j=0; j<ndetect; j++ ) {
	    data[nsample+j*MAXSAMPSCAN] = rec[inOTF->myDesc->ilocdata+j*incdatawt];
	    someOK = someOK || (data[nsample+j*MAXSAMPSCAN]!=fblank);
	  }
	}   /* end non VEGAS */
	someData = someData || someOK;  /* Any valid data? */
      } else { /* too low el - blank */
	for (j=0; j<ndetect; j++ ) data[nsample+j*MAXSAMPSCAN] = fblank;
      }
      nsample++; /* how many samples in buffer */
      
      if (lastScan<0.0) {
	lastScan   = rec[inOTF->myDesc->ilocscan]; /* Which scan number */
      }
      lastTarget = rec[inOTF->myDesc->iloctar];  /* Which target number */
      rec += inOTF->myDesc->lrec; /* Update data record pointer */
    
    } /* end loop over buffer load */
  } /* end loop reading data */
  
  /* Do last scan still in arrays*/
  if (nsample>1) {

    /* Any good data to filter? */
    if (someOK) {

      /* Filter data */
      /* Get optimum filter length - zero pad */
      nTime = ObitFFTSuggestSize (nsample*2);
      
      /* Create or resize filter as needed */
      if (filter==NULL) {
	filter = newObitTimeFilter("Cal Filter", nTime, ndetect);
      } else {
	ObitTimeFilterResize (filter, nTime);
      }
      
      /* copy data */
      for (j=0; j<ndetect; j++) {
	for (i=0; i<nsample; i++) filter->timeData[j][i] = data[j*MAXSAMPSCAN+i];
      }
      
      /* clip data */
      if (clip<1.0e19) {
	for (j=0; j<ndetect; j++) {
	  for (k=0; k<nsample; k++) {
	    if (filter->timeData[j][k]!=fblank) {
	      if (filter->timeData[j][k] > clip) filter->timeData[j][k] = clip;
	      if (filter->timeData[j][k] <-clip) filter->timeData[j][k] =-clip;
	    }
	  }
	}
      }

      /* blank fill remainder of time series */
      for (j=0; j<ndetect; j++) {
	for (i=nsample; i<nTime; i++) filter->timeData[j][i] = fblank;
      }
      
      /* Transform to frequency */
      ObitTimeFilter2Freq (filter);
      
      /* Apply filter */
      totalTime = (time[nsample-1] - time[0]);  /* Total time in scan */
      samp = totalTime / nsample;               /* sampling interval */
      parms[0] = 2.0 * samp / MAX(solInt,samp); /* Cutoff frequency in units of highest */
      ObitTimeFilterFilter (filter, -1, OBIT_TimeFilter_LowPass, parms, err);
      if (err->error) {
	Obit_log_error(err, OBIT_InfoErr, "%s: Total time %f sample time %f No. %d solInt %f", 
		       routine,totalTime, samp, nsample, solInt);
	Obit_traceback_val (err, routine, inOTF->name, outSoln);
      }
     
      /* FT back to time */
      ObitTimeFilter2Time (filter);
    } /* end of filter good data */
    
    /* Loop over solutions in scan */
    i = 0;
    t0 = time[0];
    while (i<nsample) {
      
      /* Use values in filtered time series separated by solInt */
      /* Find next time */
      while ((time[i]<t0) && (i<nsample)) {i++;}
      i = MIN (i, (nsample-1));
      
      /* Set descriptive info on Row */
      row->Time  = time[i];      /* time */
      row->TimeI = solInt;  /* time interval*/
      row->Target = (oint)(lastTarget+0.5);
      
      /* Write Soln table */
      iRow = -1;
      /* Set offset -filtered values */
      if (someOK) {  /* Some data OK */
	for (j=0; j<ndetect; j++)  row->add[j] = -filter->timeData[j][i];
      } else { /* no good data - blank */
	for (j=0; j<ndetect; j++)  row->add[j] = fblank;
      }
      if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "%s ERROR writing OTFSoln Table file", routine);
	return outSoln;
      }
      
      if (i>=(nsample-1)) break; /* done? */
      t0 += 0.25*solInt; /* time of next entry in Soln table */
      
    } /* end loop over times in scan */
  } /* end finish up data in arrays */

  
  /* Close output cal table */
  if ((ObitTableOTFSolnClose (outSoln, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s ERROR closing output OTFSoln Table file", routine);
    return outSoln;
  }
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Give warning if no data selected */
  if (!someData) Obit_log_error(err, OBIT_InfoWarn, 
				"%s: Warning: NO data selected", routine);

  /* Cleanup */
  filter = ObitTimeFilterUnref(filter);
  row = ObitTableOTFSolnUnref(row);
  if (time) g_free(time);
  if (data) g_free(data);

  return outSoln;
} /* end ObitOTFGetSolnFilter */

/**
 * Fits polynomial "baseline" to median filtered scan data as additive term in
 * an output OTFSoln table.
 * Calibration parameters are on the inOTF info member.
 * \li "solInt"    OBIT_float (1,1,1) Solution interval in days [def 10 sec].
 *                 This should not exceed 5000 samples.  Solutions will be truncated
 *                 at this limit.  This should be a sub multiple of the scan length.
 * \li "minEl"     OBIT_float (1,1,1) Minimum elevation allowed (deg)
 * \li "Order"     OBIT_int   (1,1,1) Order of polynomial to fit [def 1]
 *                 Must not exceed the number of solution intervals in a scan.
 * \param inOTF    Input OTF data. 
 * \param outOTF   OTF with which the output OTFSoln is to be associated
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created OTFSoln object which is 
 * associated with inOTF.
 */
ObitTableOTFSoln* ObitOTFGetSolnPolyBL (ObitOTF *inOTF, ObitOTF *outOTF, ObitErr *err)
{
#define MAXBLSAMPLE 1000    /* Maximum number of samples in a solution interval */
#define MAXSISAMPLE 1000    /* Maximum number of solution intervals in a scan */
  ObitTableOTFSoln *outSoln=NULL;
  ObitTableOTFSolnRow *row=NULL;
  ObitTimeFilter *filter = NULL;
  ObitOTFDesc *desc=inOTF->myDesc;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOAccess access;
  ObitInfoType type;
  ofloat lastScan=-1.0, lastTarget=-1.0;
  ofloat *rec, solInt, minEl, el, t0, tbeg=-1.0, tend=-1.0, fblank = ObitMagicF();
  ofloat *time, *data, *mtime, *mdata, *mwt, *poly;
  olong order, torder, iRow, ver, i, j, k, ibuf, lrec, nsample, msample;
  olong  npoly, ndetect, incdatawt;
  olong nfeed, nstok, nchan, nscal, ifeed, istok, incs, incfeed, ilocdata;
  ObitIOCode retCode;
  gboolean flag, doCalSelect, someOK, done, someData;
  gchar *tname;
  gchar *routine = "ObitOTFGetSolnPolyBL";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outSoln;
  g_assert (ObitOTFIsA(inOTF));

  /* Calibration wanted? */ 
  doCalSelect = FALSE;
  ObitInfoListGetTest(inOTF->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

  /* open OTF data if not already open */
  if ((inOTF->myStatus==OBIT_Inactive) || (inOTF->myStatus==OBIT_Defined)) {
    retCode = ObitOTFOpen (inOTF, access, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
  }

  /* Create output */
  tname = g_strconcat ("Calibration for: ",inOTF->name, NULL);
  ver = 0;
  npoly = 1;
  outSoln = newObitTableOTFSolnValue(tname, (ObitData*)outOTF, &ver, OBIT_IO_WriteOnly,  
				   inOTF->geom->numberDetect, npoly, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Get parameters for calibration */
  /* Solution interval default 10 sec */
  solInt = 10.0 / 86400.0;
  ObitInfoListGetTest(inOTF->info, "solInt", &type, dim, (gpointer*)&solInt);

  /* minimum allowed elevation for solution */
  minEl = 1.0;
  ObitInfoListGetTest(inOTF->info, "minEl",  &type, dim, &minEl);

  /* polynomial order */
  order = 1;
  ObitInfoListGetTest(inOTF->info, "Order",  &type, dim, &order);

  /* How many detectors? */
  if (desc->OTFType==OBIT_GBTOTF_VEGAS) {
    incs      = desc->incs;
    incfeed   = desc->incfeed;
    ilocdata  = desc->ilocdata;
    nfeed     = desc->inaxes[desc->jlocfeed];
    nstok     = desc->inaxes[desc->jlocs];
    nscal     = MIN (2, nstok);              /* Number of stokes in cal table */
    nchan     = desc->inaxes[desc->jlocf];
    ndetect   = nscal*nfeed;
    /* Better be frequency averaged */
    Obit_retval_if_fail ((nchan==1), err, outSoln,
			 "%s VEGAS data must be frequency averaged",  
			 routine);  
  } else { /* Non VEGAS */
    ndetect = inOTF->geom->numberDetect;  /* number of detectors */
  }
 /* Create work arrays */
  lrec = inOTF->myDesc->lrec;
  time  = g_malloc0(MAXBLSAMPLE*sizeof(ofloat));
  data  = g_malloc0(ndetect*MAXBLSAMPLE*sizeof(ofloat));
  mtime = g_malloc0(MAXSISAMPLE*sizeof(ofloat));
  mdata = g_malloc0(ndetect*MAXSISAMPLE*sizeof(ofloat));
  mwt   = g_malloc0(MAXSISAMPLE*sizeof(ofloat));
  poly  = g_malloc0(ndetect*(order+1)*sizeof(ofloat));
  nsample = 0;
  msample = 0;
  t0   = -1.0e20;
  lastScan = -1000.0;
  incdatawt = inOTF->myDesc->incdatawt; /* increment in data-wt axis */

  /* Open output table */
  if ((ObitTableOTFSolnOpen (outSoln, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "%s ERROR opening input OTFSoln table", routine);
    return outSoln;
  }

  /* Create Row */
  row = newObitTableOTFSolnRow (outSoln);

  /* Attach row to output buffer */
  ObitTableOTFSolnSetRow (outSoln, row, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Initialize solution row */
  row->dAz    = 0.0;     /* no pointing corrections */
  row->dEl   = 0.0;     /* no pointing corrections */
  row->TimeI  = solInt;  /* time interval*/
  row->Target = 0;
  for (j=0; j<ndetect; j++) row->mult[j] = 1.0;
  for (j=0; j<ndetect; j++) row->wt[j]   = 1.0;
  for (j=0; j<ndetect; j++) row->cal[j]  = 0.0;
  for (j=0; j<ndetect; j++) row->add[j]  = 0.0;
  for (i=0; i<npoly; i++)   row->poly[i] = 0.0;
  flag     = FALSE;
  someOK   = FALSE;
  someData = FALSE;

  /* loop over input data */
  retCode = OBIT_IO_OK;
  while (retCode == OBIT_IO_OK) {

    /* read buffer */
    if (doCalSelect) retCode = ObitOTFReadSelect (inOTF, NULL, err);
    else retCode = ObitOTFRead (inOTF, NULL, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
    done = retCode==OBIT_IO_EOF;
    if (done) break; /* done? */

    /* Record pointer */
    rec = inOTF->buffer;
  
    if (lastScan<0.0) {
      lastScan   = rec[inOTF->myDesc->ilocscan]; /* Which scan number */
      lastTarget = rec[inOTF->myDesc->iloctar];  /* Which target number */
      t0 = rec[inOTF->myDesc->iloct];            /* Start time */
    }

    /* Loop over buffer */
    for (ibuf=0; ibuf<inOTF->myDesc->numRecBuff; ibuf++) {

      /* Finished solution interval? */
      if ((rec[inOTF->myDesc->iloct] > (t0+solInt)) || 
	  /* or end of scan */
	  ( rec[inOTF->myDesc->ilocscan] != lastScan) ||
	  /* or exceed limit on number of samples */
	  (nsample>=MAXBLSAMPLE)) {
	if (msample==0) tbeg = t0;  /* save start time */
	tend = time[nsample-1];     /* save end time */
	mtime[msample] = meanValue(time, 1, nsample) - tbeg;
	mwt[msample] = (ofloat)nsample;   /* weight by number of samples */
	if (someOK) {
	  for (j=0; j<ndetect; j++ )
	    mdata[msample+j*MAXSISAMPLE] = 
	      GetSolnMedianAvg(&data[j*MAXBLSAMPLE], 1, nsample);
	} else {
	  for (j=0; j<ndetect; j++ ) mdata[msample+j*MAXSISAMPLE] = fblank;
	}
	if (msample<MAXSISAMPLE-1) msample++;
	/* reset SI accumulators */
	nsample = 0;
	t0 = rec[inOTF->myDesc->iloct];
	someOK = FALSE;
      } /* end process SI */

      /* Finished Scan? */
      if (rec[inOTF->myDesc->ilocscan] != lastScan) {
	/* Fit polynomials */
	if (msample>order) torder = order;
	else torder = MAX (0, msample-1);
	for (j=0; j<ndetect; j++ )
	  FitBLPoly (&poly[j*(torder+1)], torder, mtime, &mdata[j*MAXSISAMPLE],
		     mwt, msample);

	/* Write beginnning (tbeg) solution */
	/* Fill Soln row */
	row->Target = (oint)(lastTarget+0.5);
	row->Time  = tbeg;
	for (j=0; j<ndetect; j++)
	  row->add[j] = -EvalPoly (torder, &poly[j*(torder+1)], tbeg - tbeg);
	iRow = -1;
	ObitTableOTFSolnWriteRow (outSoln, iRow, row, err);
	if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

	/* Write output solutions */
	for (k=0; k<msample; k++) {
	  /* Fill Soln row */
	  row->Time  = mtime[k] + tbeg;
	  for (j=0; j<ndetect; j++)
	    row->add[j] = -EvalPoly (torder, &poly[j*(torder+1)],  mtime[k]);
	  ObitTableOTFSolnWriteRow (outSoln, iRow, row, err);
	  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
	} /* end loop over solution intervals */

	/* Add end (tend) solution */
	/* Fill Soln row */
	row->Time  = tend;
	for (j=0; j<ndetect; j++)
       	  row->add[j] = -EvalPoly (torder, &poly[j*(torder+1)], tend - tbeg); 
	ObitTableOTFSolnWriteRow (outSoln, iRow, row, err);
	if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

	/* reset scan accumulators */
	lastScan = rec[inOTF->myDesc->ilocscan];
	msample = 0;
      } /* end finish scan */

      /* Is data below the elevation limit? */
      el = ObitOTFArrayGeomElev (inOTF->geom, rec[inOTF->myDesc->iloct], 
				 rec[inOTF->myDesc->ilocra], 
				 rec[inOTF->myDesc->ilocdec]);
      flag = (el < minEl);
      
      /* accumulate */
      time[nsample] = rec[inOTF->myDesc->iloct];
      if (!flag) { /* OK */
	/* VEGAS different use, parallel poln and feeds */
	j = 0;
	if (desc->OTFType==OBIT_GBTOTF_VEGAS) {
	  /* Loop over feed */
	  for (ifeed=0; ifeed<nfeed; ifeed++) {
	    /* Loop over parallel Stokes */
	    for (istok=0; istok<nscal; istok++) {
	       data[nsample+j*MAXSAMPLE] = 
		 rec[ilocdata + ifeed*incfeed + istok*incs];
	       someOK = someOK || (data[nsample+j*MAXBLSAMPLE]!=fblank);
	       j++;
	    } /* end Stokes loop */
	  } /* end feed loop */
	} else { /* Non VEGAS */
	  for (j=0; j<ndetect; j++ ) {
	    data[nsample+j*MAXBLSAMPLE] = rec[inOTF->myDesc->ilocdata+j*incdatawt];
	    someOK = someOK || (data[nsample+j*MAXBLSAMPLE]!=fblank);
	  }
	}   /* end non VEGAS */
	someData = someData || someOK;  /* Any valid data? */
      } else { /* too low el - blank */
	for (j=0; j<ndetect; j++ ) data[nsample+j*MAXBLSAMPLE] = fblank;
      }
      nsample++; /* how many samples in buffer */
       
      lastTarget = rec[inOTF->myDesc->iloctar];  /* Which target number */
      rec += inOTF->myDesc->lrec; /* Update data record pointer */
      
    } /* end loop over buffer load */
  } /* end loop reading data */
  
  /* Do last SI still in arrays*/
  if (nsample>1) {
    /* Finished solution interval? */
    if (msample==0) tbeg = t0;  /* save start time */
    tend = time[nsample-1];     /* save end time */
    mtime[msample] = meanValue(time, 1, nsample) - tbeg;
    mwt[msample] = (ofloat)nsample;   /* weight by number of samples */
    if (someOK) {
      for (j=0; j<ndetect; j++ ) 
	mdata[msample+j*MAXSISAMPLE] = GetSolnMedianAvg(&data[j*MAXBLSAMPLE], 1, nsample);
    } else {
      for (j=0; j<ndetect; j++ ) mdata[msample+j*MAXSISAMPLE] = fblank;
    }
    if (msample<MAXSISAMPLE-1) msample++;
  } /* end process SI */
    
  /* Something left in last scan? */
  if (msample>1) {
    /* Fit polynomials */
    if (msample>order) torder = order;
    else torder = MAX (0, msample-1);
    for (j=0; j<ndetect; j++ ) 
      FitBLPoly (&poly[j*(torder+1)], torder, mtime, &mdata[j*MAXSISAMPLE], 
		 mwt, msample); 
    
    /* Write beginnning (tbeg) solution */
    /* Fill Soln row */
    row->Target = (oint)(lastTarget+0.5);
    row->Time  = tbeg;
    for (j=0; j<ndetect; j++)
      row->add[j] = -EvalPoly (torder, &poly[j*(torder+1)], tbeg - tbeg);
    iRow = -1;
    ObitTableOTFSolnWriteRow (outSoln, iRow, row, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
    
    /* Write output solutions */
    for (k=0; k<msample; k++) {
      /* Fill Soln row */
      row->Time  = mtime[k] + tbeg;
      for (j=0; j<ndetect; j++) 
	row->add[j] = -EvalPoly (torder, &poly[j*(torder+1)],  mtime[k]);
      ObitTableOTFSolnWriteRow (outSoln, iRow, row, err);
      if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
    } /* end loop over solution intervals */
    
    /* Add end (tend) solution */
    /* Fill Soln row */
    row->Time  = tend;
    for (j=0; j<ndetect; j++) 
      row->add[j] = -EvalPoly (torder, &poly[j*(torder+1)], tend - tbeg);
    ObitTableOTFSolnWriteRow (outSoln, iRow, row, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
  } /* end finish scan */

  
  /* Close output cal table */
  ObitTableOTFSolnClose (outSoln, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Give warning if no data selected */
  if (!someData) Obit_log_error(err, OBIT_InfoWarn, 
				"%s: Warning: NO data selected", routine);

  /* Cleanup */
  filter = ObitTimeFilterUnref(filter);
  row = ObitTableOTFSolnUnref(row);
  if (time)  g_free(time);
  if (data)  g_free(data);
  if (mtime) g_free(mtime);
  if (mdata) g_free(mdata);
  if (mwt)   g_free(mwt);
  if (poly ) g_free(poly);

  return outSoln;
} /* end ObitOTFGetSolnPolyBL */

/**
 * Baseline removal for multibeam instrument.
 * Fit one term, time variable common, atmospheric polynomial and a single offset
 * per detector.
 * Since the different detectors each have an individual multiplicative term, the 
 * Atmospheric + offset are places in the the detector's additive term and the
 * polynomical is set to zero.
 * Scans in excess of 5000 samples will be broken into several.
 * Either detector offsets, common mode variations or both may be selected to Soln.
 * Calibration parameters are on the inOTF info member.
 * \li "solInt"    OBIT_float (1,1,1) Solution interval in days [def 10 sec].
 * \li "maxInt"    OBIT_float (1,1,1) max. Interval in days [def 10 min].
 *                 Scans longer than this will be broken into pieces
 * \li "Tau0"      OBIT_float (?,?,1) Zenith opacity in nepers [def 0].
 *                 Can be passed as either a constant scalar or
 *                 an array of (time, tau0) pairs.
 * \li "minEl"     OBIT_float (1,1,1) Minimum elevation allowed (deg)
 * \li "clipSig"   OBIT_float (1,1,1) data outside of range +/- Clip sigma are blanked
 *                                    [def very large -> no clipping]
 * \li "plotDet"   OBIT_long   (1,1,1) Detector number to plot per scan [def =-1 = none]
 * \li "calType"   OBIT_string (*,1,1) Calibration type desired [def "Both"]
 *                 "Both" => Detector offsets and common mode poly
 *                 "Common" => common mode poly
 *                 "Offset" Detector offsets 
 * \param inOTF    Input OTF data. 
 * \param outOTF   OTF with which the output OTFSoln is to be associated
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created OTFSoln object which is 
 * associated with inOTF.
 */
ObitTableOTFSoln* ObitOTFGetSolnMBBase (ObitOTF *inOTF, ObitOTF *outOTF, ObitErr *err)
{
  ObitTableOTFSoln *outSoln=NULL;
  ObitTableOTFSolnRow *row=NULL;
  gint32 taudim[MAXINFOELEMDIM] = {1,1,1,1,1}, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ObitIOAccess access;
  ofloat *tau0, tMult, airmass;
  ofloat lastScan=-1.0, lastTarget=-1.0, lastTime=-1000.0;
  ofloat *rec=NULL, solInt, maxInt, minEl, el, t0, fblank = ObitMagicF();
  ofloat *time=NULL, *data=NULL, *poly=NULL, *tpoly=NULL, *offset=NULL, *wt=NULL, clipsig;
  olong iRow, ver, i, j, k, ibuf, lrec, nsample, npoly, off;
  olong  mpoly, ndetect, plotDetect, incdatawt;
  ObitIOCode retCode;
  gboolean flag, doCalSelect, someOK, someData, done, doCommon, doOffset;
  gchar *tname, calType[50];
  gchar *routine = "ObitOTFGetSolnMBBase";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outSoln;
  g_assert (ObitOTFIsA(inOTF));

   /* Calibration wanted? */ 
  doCalSelect = FALSE;
  ObitInfoListGetTest(inOTF->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

  /* open OTF data if not already open */
  if ((inOTF->myStatus==OBIT_Inactive) || (inOTF->myStatus==OBIT_Defined)) {
    retCode = ObitOTFOpen (inOTF, access, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
  }
  
  /* Create output */
  tname = g_strconcat ("Calibration for: ",inOTF->name, NULL);
  ver = 0;
  mpoly = 1;
  outSoln = newObitTableOTFSolnValue(tname, (ObitData*)outOTF, &ver, OBIT_IO_WriteOnly,  
				     inOTF->geom->numberDetect, mpoly, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Create work arrays */
  ndetect = inOTF->geom->numberDetect;  /* number of detectors */
  lrec    = inOTF->myDesc->lrec;
  time    = g_malloc0(MAXSAMPSCAN*sizeof(ofloat));
  data    = g_malloc0(ndetect*MAXSAMPSCAN*sizeof(ofloat));
  wt      = g_malloc0(ndetect*MAXSAMPSCAN*sizeof(ofloat));
  offset  = g_malloc0(ndetect*sizeof(ofloat));
  nsample = 0;
  t0      = -1.0e20;  lastScan = -1000.0;
  time[0] = 1.0e20;                     /* No data yet */
  incdatawt = inOTF->myDesc->incdatawt; /* increment in data-wt axis */

  /* Get parameters for calibration */
  /* Solution interval default 10 sec */
  solInt = 10.0 / 86400.0;
  ObitInfoListGetTest(inOTF->info, "solInt", &type, dim, (gpointer*)&solInt);

  /* maximum interval default 10 min */
  maxInt = 10.0 * 60.0 / 86400.0;
  ObitInfoListGetTest(inOTF->info, "maxInt", &type, dim, (gpointer*)&maxInt);

  /* minimum allowed elevation for solution */
  minEl = 1.0;
  ObitInfoListGetTest(inOTF->info, "minEl",  &type, dim, (gpointer*)&minEl);

  /* Clip value */
  clipsig = 1.0e20;
  ObitInfoListGetTest(inOTF->info, "clipSig",  &type, dim, (gpointer*)&clipsig);
  if (clipsig<1.0e19) 
    Obit_log_error(err, OBIT_InfoErr, "%s: Clipping outliers at %7.2f sigma", 
		   routine, clipsig);

  /* Opacity */
  tau0 = NULL;
  ObitInfoListGetP(inOTF->info, "Tau0",   &type, taudim, (gpointer*)&tau0);

  /* Detector to plot */
  plotDetect = -1;
  ObitInfoListGetTest(inOTF->info, "plotDet",  &type, dim, &plotDetect);

  /* calibration type */
  for (i=0; i<50; i++) calType[i] = 0;
  strcpy (calType, "Both");
  ObitInfoListGetTest(inOTF->info, "calType", &type, dim, calType);
  /* What calibrations are desired? */
  doCommon = FALSE; doOffset = FALSE;
  if (!strncmp(calType, "Common",4)) doCommon = TRUE;
  else if (!strncmp(calType, "Offset",6)) doOffset = TRUE;
  else if (!strncmp(calType, "Both",10)) {doCommon=TRUE; doOffset = TRUE;}
  else {doCommon=TRUE; doOffset = TRUE;}

  /* Open output table */
  if ((ObitTableOTFSolnOpen (outSoln, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "%s ERROR opening input OTFSoln table", routine);
    goto cleanup;
  }

  /* Create Row */
  row = newObitTableOTFSolnRow (outSoln);

  /* Attach row to output buffer */
  ObitTableOTFSolnSetRow (outSoln, row, err);
  if (err->error) goto cleanup;

  /* Initialize solution row */
  row->dAz = 0.0;  /* no pointing corrections */
  row->dEl = 0.0;  /* no pointing corrections */
  row->Target = 0;
  for (j=0; j<ndetect; j++) row->mult[j] = 1.0;
  for (j=0; j<ndetect; j++) row->wt[j]   = 1.0;
  for (j=0; j<ndetect; j++) row->cal[j]  = 0.0;
  for (j=0; j<ndetect; j++) row->add[j]  = 0.0;
  flag     = FALSE;
  someOK   = FALSE;
  someData = FALSE;
  
  /* loop over input data */
  retCode = OBIT_IO_OK;
  while (retCode == OBIT_IO_OK) {

    /* read buffer */
    if (doCalSelect) retCode = ObitOTFReadSelect (inOTF, NULL, err);
    else retCode = ObitOTFRead (inOTF, NULL, err);
    if (err->error) goto cleanup;
    done = retCode==OBIT_IO_EOF;
    if (done) break; /* done? */

    /* Record pointer */
    rec = inOTF->buffer;
  
    if (lastScan<0.0) {
      lastScan   = rec[inOTF->myDesc->ilocscan]; /* Which scan number */
      lastTarget = rec[inOTF->myDesc->iloctar];  /* Which target number */
    }

    /* Loop over buffer */
    for (ibuf=0; ibuf<inOTF->myDesc->numRecBuff; ibuf++) {

      /* Scan read? maxInt? If so, compute solution and write. Also if work arrays full */
      if (( rec[inOTF->myDesc->ilocscan] != lastScan)  || /* new scan */
	  ((rec[inOTF->myDesc->iloct]-time[0])>maxInt) || /* or maxInt Done */
	  (nsample>=MAXSAMPSCAN) || done) {   /* or exceed limit on number of samples, or done */
	
	/* Any good data to process? */
	if (someOK) {
	  
	  /* time relative to first integration */
	  t0 = time[0];
	  for (i=0; i<nsample; i++) time[i] -= t0;
      
	  /* DEBUG
	  fprintf (stderr,"DEBUG lastScan %f \n", lastScan); */

	  /* Fit data */
	  FitMBBLPoly (solInt, &npoly, &tpoly, &poly, offset, ndetect, MAXSAMPSCAN, time, 
		       data, wt, nsample);

	  /* reject outlyers */
	  FitMBBLOut (npoly, tpoly, poly, offset, ndetect, MAXSAMPSCAN, time, 
		      data, wt, nsample, clipsig);
	  /* Again to be sure */
	  FitMBBLOut (npoly, tpoly, poly, offset, ndetect, MAXSAMPSCAN, time, 
		      data, wt, nsample, clipsig);
	  if (tpoly) g_free(tpoly); tpoly = NULL;/* free work memory */
	  if (poly)  g_free(poly);  poly  = NULL;

	  /* refit data */
	  FitMBBLPoly (solInt, &npoly, &tpoly, &poly, offset, ndetect, MAXSAMPSCAN, time, 
		       data, wt, nsample);

	  /* Plot? */
	  if (plotDetect>0) 
	    PlotMBBL (npoly, tpoly, poly, offset, ndetect, MAXSAMPSCAN, time, 
		      data, wt, nsample, plotDetect, t0, err);
	  if (err->error) goto cleanup;
	} /* end of fit good data */
	
	/* Loop over time segments */
	off = 0;
	if (tpoly==NULL) npoly = 0;  /* In case no good data */
	for (k=0; k<npoly; k++) {

	  /* Set descriptive info on Row */
	  row->TimeI = solInt;  /* time interval*/
	  row->Target = (oint)(lastTarget+0.5);
	  
	  row->Time  = tpoly[k] + t0;      /* start time */

	  /* Opacity correction */
	  el = ObitOTFArrayGeomElev (inOTF->geom, row->Time, 
				     rec[inOTF->myDesc->ilocra], 
				     rec[inOTF->myDesc->ilocdec]);
	  if (el > 0.1) 
	    airmass = 1.0 / cos (DG2RAD * (90.0 - el));
	  else
	    airmass = 10.0;
	  tMult = exp (getTau0(row->Time,taudim,tau0) * airmass); /* multiplicative term */
	  for (j=0; j<ndetect; j++) row->mult[j] = tMult;
	
	  /* solution (as correction) */
	  if (someOK && (poly[2*k]!=fblank)) {  /* Some data and fit OK */
	    /* Find and evaluate time segment */
	    while ((tpoly[k+1]>time[off]) && (off<nsample)) {off++;}
	    if (doCommon) row->poly[0] = -(poly[2*k] + poly[2*k+1]*(tpoly[k]));
	    else row->poly[0] = 0.0;
	    /* row->poly[0] = 0.0;  DEBUG */
	    for (j=0; j<ndetect; j++)  {
	      if (doOffset) {
		if (offset[j]!=fblank) row->add[j] = -offset[j];
		else row->add[j] = fblank;
	      } else  row->add[j] = 0.0;
	    }
	  } else { /* no good data - blank */
	    for (j=0; j<ndetect; j++)  row->add[j] = fblank;
	    row->poly[0] = fblank;
	  }

	  /* Write Soln table for beginning and end of segment */
	  iRow = -1;
	  if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	       != OBIT_IO_OK) || (err->error>0)) { 
	    Obit_log_error(err, OBIT_Error, "%s ERROR writing OTFSoln Table file", routine);
	    goto cleanup;
	  }
	  row->Time  = MIN (t0, lastTime);      /* end time */

	  /* Opacity correction */
	  el = ObitOTFArrayGeomElev (inOTF->geom, row->Time, 
				     rec[inOTF->myDesc->ilocra], 
				     rec[inOTF->myDesc->ilocdec]);
	  if (el > 0.1) 
	    airmass = 1.0 / cos (DG2RAD * (90.0 - el));
	  else
	    airmass = 10.0;
	  tMult = exp (getTau0(row->Time,taudim,tau0) * airmass); /* multiplicative term */
	  for (j=0; j<ndetect; j++) row->mult[j] = tMult;
	
	  /* solution (as correction) */
	  if (someOK && (poly[2*k]!=fblank)) {  /* Some data and fit OK */
	    row->Time  = MIN (tpoly[k+1]+t0, lastTime);      /* end time */
	    /* Find and evaluate time segment */
	    while ((tpoly[k+1]>time[off]) && (off<nsample)) {off++;}
	    row->poly[0] = -(poly[2*k] + poly[2*k+1]*(tpoly[k+1]));
	    for (j=0; j<ndetect; j++)  {
	      if (offset[j]!=fblank) row->add[j] = -offset[j];
	      else row->add[j] = fblank;
	    }
	  } else { /* no good data - blank */
	    for (j=0; j<ndetect; j++)  row->add[j] = fblank;
	    row->poly[0] = fblank;
	  }

	  if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	       != OBIT_IO_OK) || (err->error>0)) { 
	    Obit_log_error(err, OBIT_Error, "%s ERROR writing OTFSoln Table file", routine);
	    goto cleanup;
	  }
	} /* end loop over time segments */
	
	/* initialize accumulators */
	if (tpoly) g_free(tpoly);  tpoly = NULL;/* free work memory */
	if (poly)  g_free(poly);   poly  = NULL;
	nsample = 0;
	lastScan = -1000.0;
	someOK = FALSE;
	
      } /* end do solution */
      
      /* Is data below the elevation limit? */
      el = ObitOTFArrayGeomElev (inOTF->geom, rec[inOTF->myDesc->iloct], 
				 rec[inOTF->myDesc->ilocra], 
				 rec[inOTF->myDesc->ilocdec]);
      flag = (el < minEl);

      /* accumulate */
      time[nsample] = rec[inOTF->myDesc->iloct];
      lastTime = time[nsample];
      if (!flag) { /* OK */
	for (j=0; j<ndetect; j++ ) {
	  data[nsample+j*MAXSAMPSCAN] = rec[inOTF->myDesc->ilocdata+j*incdatawt];
	  someOK = someOK || (data[nsample+j*MAXSAMPSCAN]!=fblank);
	  if (data[nsample+j*MAXSAMPSCAN]!=fblank) {
	    wt[nsample+j*MAXSAMPSCAN] = 1.0;
	  } else {
	    wt[nsample+j*MAXSAMPSCAN] = 0.0;
	  }
	}
	someData = someData || someOK;  /* Any valid data? */
      } else { /* too low el - blank */
	for (j=0; j<ndetect; j++ ) data[nsample+j*MAXSAMPSCAN] = fblank;
	for (j=0; j<ndetect; j++ ) wt[nsample+j*MAXSAMPSCAN] = 0.0;
      }
      nsample++; /* how many samples in buffer */
      
      if (lastScan<0.0) {
	lastScan   = rec[inOTF->myDesc->ilocscan]; /* Which scan number */
      }
      lastTarget = rec[inOTF->myDesc->iloctar];  /* Which target number */
      rec += inOTF->myDesc->lrec; /* Update data record pointer */
      
    } /* end loop over buffer load */
  } /* end loop reading data */
  
  /* Do last scan still in arrays*/
  if (nsample>1) {

    /* Any good data to process? */
    if (someOK) {

      /* time relative to first integration */
      t0 = time[0];
      for (i=0; i<nsample; i++) time[i] -= t0;
      
      /* DEBUG
      fprintf (stderr,"DEBUG lastScan %f \n", lastScan); */
      /* Fit data */
      FitMBBLPoly (solInt, &npoly, &tpoly, &poly, offset, ndetect, MAXSAMPSCAN, time, 
		   data, wt, nsample);
      
      /* reject outlyers */
      FitMBBLOut (npoly, tpoly, poly, offset, ndetect, MAXSAMPSCAN, time, 
		  data, wt, nsample, clipsig);
      /* Again to be sure */
      FitMBBLOut (npoly, tpoly, poly, offset, ndetect, MAXSAMPSCAN, time, 
		  data, wt, nsample, clipsig);
      if (tpoly) g_free(tpoly);  tpoly = NULL;/* free work memory */
      if (poly)  g_free(poly);    poly = NULL;
      
      /* refit data */
      FitMBBLPoly (solInt, &npoly, &tpoly, &poly, offset, ndetect, MAXSAMPSCAN, time, 
		   data, wt, nsample);
      
      /* Plot? */
      if (plotDetect>0) 
	PlotMBBL (npoly, tpoly, poly, offset, ndetect, MAXSAMPSCAN, time, 
		  data, wt, nsample, plotDetect, t0, err);
      if (err->error) goto cleanup;
    } /* end of fit good data */
    
    /* Loop over time segments */
    off = 0;
    if (tpoly==NULL) npoly = 0;  /* In case no good data */
    for (k=0; k<npoly; k++) {
      
      /* Set descriptive info on Row */
      row->TimeI = 0.5*solInt;     /* time interval*/
      row->Target = (oint)(lastTarget+0.5);
      
      row->Time  = tpoly[k] + t0;      /* start time */
      /* solution (as correction) */
      if (someOK && (poly[2*k]!=fblank)) {  /* Some data OK */
	/* Find and evaluate time segment */
	while ((tpoly[k+1]>time[off]) && (off<nsample)) {off++;}
	if (doCommon) row->poly[0] = -(poly[2*k] + poly[2*k+1]*(tpoly[k]));
	else row->poly[0] = 0.0;
	for (j=0; j<ndetect; j++)  {
	  if (doOffset) {
	    if (offset[j]!=fblank) row->add[j] = -offset[j];
	    else row->add[j] = fblank;
	  } else  row->add[j] = 0.0;
	}
      } else { /* no good data - blank */
	for (j=0; j<ndetect; j++)  row->add[j] = fblank;
	row->poly[0] = fblank;
      }
      
      /* Write Soln table for beginning and end of segment */
      iRow = -1;
      if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "%s ERROR writing OTFSoln Table file", routine);
	goto cleanup;
      }
      row->Time  = MIN (t0, lastTime);      /* end time */

      /* Opacity correction */
      el = ObitOTFArrayGeomElev (inOTF->geom, row->Time, 
				 rec[inOTF->myDesc->ilocra], 
				 rec[inOTF->myDesc->ilocdec]);
      if (el > 0.1) 
	airmass = 1.0 / cos (DG2RAD * (90.0 - el));
      else
	airmass = 10.0;
      tMult = exp (getTau0(row->Time,taudim,tau0) * airmass); /* multiplicative term */
      for (j=0; j<ndetect; j++) row->mult[j] = tMult;
	
      /* solution (as correction) */
      if (someOK && (poly[2*k]!=fblank)) {  /* Some data OK */
	row->Time  = MIN (tpoly[k+1]+t0, lastTime);      /* end time */
	/* Find and evaluate time segment */
	while ((tpoly[k+1]>time[off]) && (off<nsample)) {off++;}
	row->poly[0] = -(poly[2*k] + poly[2*k+1]*(tpoly[k+1]));
	/* row->poly[0] = 0.0;  DEBUG */
	for (j=0; j<ndetect; j++)  {
	  if (offset[j]!=fblank) row->add[j] = -offset[j];
	  else row->add[j] = fblank;
	}
      } else { /* no good data - blank */
	for (j=0; j<ndetect; j++)  row->add[j] = fblank;
	row->poly[0] = fblank;
      }
      if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "%s ERROR writing OTFSoln Table file", routine);
	goto cleanup;
      } /* end loop over time segments */
    } /* end loop over polynomials */
  } /* end finish up data in arrays */
    
  
  /* Close output cal table */
  if ((ObitTableOTFSolnClose (outSoln, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s ERROR closing output OTFSoln Table file", routine);
    goto cleanup;
  }
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) goto cleanup;

  /* Give warning if no data selected */
  if (!someData) Obit_log_error(err, OBIT_InfoWarn, 
				"%s: Warning: NO data selected", routine);

  /* Deallocate arrays */
  cleanup:
  row = ObitTableOTFSolnUnref(row);
  if (time)   g_free(time);
  if (data)   g_free(data);
  if (wt)     g_free(wt);
  if (offset) g_free(offset);
  if (tpoly)  g_free(tpoly);  /* free work memory */
  if (poly)   g_free(poly);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  return outSoln;
} /* end ObitOTFGetSolnMBBase */

/**
 * Determine estimates of the cal value from a set of data and then
 * calculate the median value.  
 * Successive off/on/off triplets are used for each estimate.
 * \param n      number of data values
 * \param isCal  For each datum, TRUE if cal on.
 * \param data   On input the data values, will be replaced on output
 *               by the sequence of valid estimates of the cal.
 * \return       Median Cal value or fblank if not possible
 */
ofloat ObitOTFGetSolnAvgCal (olong n, gboolean *isCal, ofloat *data)
{
  olong i, j, prior, follow, ngood=0;
  ofloat fblank = ObitMagicF();

  /* loop over data looking for cal on with previous and following cal off*/
  ngood = 0;
  for (i=1; i<n-1; i++) {

    /* if cal off ignore this one of invalid datum */
    if ((!isCal[i]) || (data[i]==fblank)) continue;

    /* find previous good cal off */
    prior = -1;
    for (j=i-1; j>=0; j--) {
      if (!isCal[j] && (data[j]!=fblank)) {
	prior = j;
	break;
      }
    }

    /* if not found, ignore */
    if (prior<0) continue;

    /* find following good cal off */
    follow = -1;
    for (j=i+1; j<n; j++) {
      if (!isCal[j] && (data[j]!=fblank)) {
	follow = j;
	break;
      }
    }

    /* if not found, ignore */
    if (follow<0) continue;

    /* Add it to the good list */
    data[ngood++] = data[i] - 0.5*(data[prior] + data[follow]);

  } /* end loop over data */

  /* any good data? */
  if (ngood<1) return fblank;

  /* Get median */
  return GetSolnMedianAvg (data, 1, ngood);
} /* end  ObitOTFGetSolnAvgCal */

/**
 * Determine instrumental calibration for an OTF multibeam OTF dataset.
 * Calculates average gain from cal measurements
 * calculates median ofset
 * Calibration parameters are on the inOTF info member.
 * \li "solInt"   OBIT_float (1,1,1) Solution interval in days [def 1 sec].
 * \param inOTF    Input OTF data. 
 * \param outOTF   OTF with which the output  OTFSoln is to be associated
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created OTFSoln object which is 
 * associated with outOTF.
 */
ObitTableOTFSoln* ObitOTFGetInstCal (ObitOTF *inOTF, ObitOTF *outOTF, ObitErr *err)
{
#define MAXSAMPLE 1000   /* Maximum number of samples in an integration */
  ObitTableOTFSoln *outSoln=NULL;
  ObitTableOTFSolnRow *row=NULL;
  ObitOTFDesc *desc=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat *rec, solInt, t0, *iCal, val, **accum, sumTime, fblank = ObitMagicF();
  ofloat lastTarget=-1.0, lastScan=-1.0, lastTime=-1.0, *lastCal=NULL;
  olong iRow, ver, i, j, k, lrec, ncoef;
  olong  npoly, nDet, nTime, incdatawt;
  ObitIOCode retCode;
  gboolean someData;
  gchar *tname;
  gchar *routine = "ObitOTFGetInstCal";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outSoln;
  g_assert (ObitOTFIsA(inOTF));
  desc = inOTF->myDesc;

  /* open OTF data to fully instantiate if not already open */
  if ((inOTF->myStatus==OBIT_Inactive) || (inOTF->myStatus==OBIT_Defined)) {
    retCode = ObitOTFOpen (inOTF, OBIT_IO_ReadWrite, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
  }

  /* Create output */
  tname = g_strconcat ("Calibration for: ",inOTF->name, NULL);
  ver = 0;
  ncoef = 1;
  npoly = ncoef;
  outSoln = newObitTableOTFSolnValue(tname, (ObitData*)outOTF, &ver, OBIT_IO_WriteOnly,  
				   inOTF->geom->numberDetect, npoly, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Create work arrays */
  nDet = inOTF->geom->numberDetect;
  lrec = inOTF->myDesc->lrec;
  iCal = g_malloc0(MAXSAMPLE*sizeof(ofloat)); 
  lastCal = g_malloc0(nDet*sizeof(ofloat)); 
  for (i=0; i<nDet; i++) lastCal[i] = fblank;
  accum = g_malloc0(nDet*sizeof(ofloat*)); 
  for (i=0; i<nDet; i++) accum[i] = g_malloc0(MAXSAMPLE*sizeof(ofloat));
  t0 = -1.0e20;
  incdatawt = inOTF->myDesc->incdatawt; /* increment in data-wt axis */
  someData = FALSE;

  /* Get parameters for calibration */
  /* Solution interval default 1 sec */
  solInt = 1.0 / 86400.0;
  ObitInfoListGetTest(inOTF->info, "solInt", &type, dim, (gpointer*)&solInt);

  /* Open table */
  if ((ObitTableOTFSolnOpen (outSoln, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input OTFSoln table");
    return outSoln;
  }

  /* Create Row */
  row = newObitTableOTFSolnRow (outSoln);

  /* Attach row to output buffer */
  ObitTableOTFSolnSetRow (outSoln, row, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* loop calibrating data */
  retCode = OBIT_IO_OK;
  sumTime = 0.0;
  nTime = 0;
  while (retCode == OBIT_IO_OK) {
    
    /* read buffer */
    retCode = ObitOTFRead (inOTF, NULL, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
    if (retCode==OBIT_IO_EOF) break; /* done? */
    
    /* Record pointer */
    rec = inOTF->buffer;
    
    /* First time */
    if (t0<-1.0e10) {
      t0 = rec[inOTF->myDesc->iloct];
      lastTime   = rec[inOTF->myDesc->iloct];
      lastTarget = rec[inOTF->myDesc->iloctar];
      lastScan   = rec[inOTF->myDesc->ilocscan];
    }
    
    /* Loop over buffer */
    for (i=0; i<inOTF->myDesc->numRecBuff; i++) {
      
      /* Accumulation finished? If so, compute solution and write.*/
      if ((rec[inOTF->myDesc->iloct] > (t0+solInt)) || (nTime>=MAXSAMPLE) ||
	  (rec[inOTF->myDesc->iloctar] != lastTarget) ||  
	  (rec[inOTF->myDesc->ilocscan] != lastScan)) {
	
	/* Not first time - assume first descriptive parameter never blanked */
	if (nTime>0) {
	  /* Set descriptive info on Row */
	  row->Time  = sumTime/nTime;  /* center time */
	  row->TimeI = 2.0 * (row->Time - t0);
	  
	  /* Get calibration */
	  ObitOTFGetInstSolve (inOTF->geom, inOTF->myDesc, nDet, nTime, accum, iCal, row);

	  /* Propagate last good cals if needed */
	  for (j=0; j<nDet; j++ ) {
	    if (row->cal[j]==fblank) {
	      row->cal[j]  = lastCal[j];
	      row->mult[j] = 1.0 / lastCal[j];
	    }
	    lastCal[j] = row->cal[j];
	  }
	  
	  /* Write Soln table - bracket interval */
	  iRow = -1;
	  row->Target = (oint)(lastTarget+0.5);
	  row->Time   = t0;
	  if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	       != OBIT_IO_OK) || (err->error>0)) { 
	    Obit_log_error(err, OBIT_Error, "ERROR writing OTFSoln Table file");
	    return outSoln;
	  }
	  row->Time   = lastTime;
	  if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	       != OBIT_IO_OK) || (err->error>0)) { 
	    Obit_log_error(err, OBIT_Error, "ERROR writing OTFSoln Table file");
	    return outSoln;
	  }
	  /* initialize accumulators */
	  t0 = rec[inOTF->myDesc->iloct];
	  lastTime   = rec[inOTF->myDesc->iloct];
	  lastTarget = rec[inOTF->myDesc->iloctar];
	  lastScan   = rec[inOTF->myDesc->ilocscan];
	  sumTime = 0.0;
	  nTime = 0;
	  for (k=0; k<MAXSAMPLE; k++) iCal[k] = 0;
	  for (j=0; j<nDet; j++ ) {
	    for (k=0; k<MAXSAMPLE; k++) accum[j][k] = 0.0;
	  }
	} /* end of do solution if there is data */
	
      } /* end do solution */
      
      /* accumulate */
      lastTime = rec[inOTF->myDesc->iloct];
      sumTime += rec[inOTF->myDesc->iloct];
      iCal[nTime] = fabs(rec[inOTF->myDesc->iloccal]);
      for (j=0; j<nDet; j++ ) {
	val = rec[desc->ilocdata+j*incdatawt];
	if (val != fblank) {
	  accum[j][nTime] = val;
	}
      }
      nTime++; /* how many data points */
      someData = TRUE;
      rec += inOTF->myDesc->lrec; /* Data record pointer */
      
    } /* end loop over buffer load */
  } /* end loop reading/gridding data */
  
  /* Finish up any data in accumulator */
  if (nTime>0) {
    /* Set descriptive info on Row */
    row->Time  = sumTime/nTime; /* time */
    row->TimeI = 2.0 * (row->Time - t0);
    
    /* Get calibration */
    ObitOTFGetInstSolve (inOTF->geom, inOTF->myDesc, nDet, nTime, accum, iCal, row);
    rec += inOTF->myDesc->lrec; /* Data record pointer */
    
    /* Propagate last good cals if needed */
    for (j=0; j<nDet; j++ ) {
      if (row->cal[j]==fblank) {
	row->cal[j]  = lastCal[j];
	row->mult[j] = 1.0 / lastCal[j];
      }
    }

    /* Write Soln table - bracket interval */
    iRow = -1;
    row->Time = t0;
    if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR writing OTFSoln Table file");
      return outSoln;
    }
    row->Time   = lastTime;
    if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR writing OTFSoln Table file");
      return outSoln;
    }
  } /* End final solution */

  /* Close cal table */
  if ((ObitTableOTFSolnClose (outSoln, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input OTFSoln Table file");
    return outSoln;
  }
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Give warning if no data selected */
  if (!someData) Obit_log_error(err, OBIT_InfoWarn, 
				"%s: Warning: NO data selected", routine);

  /* Cleanup */
  for (k=0; k<nDet; k++) if (accum[k]) g_free(accum[k]);
  if (accum)   g_free(accum);
  if (iCal)    g_free(iCal);
  if (lastCal) g_free(lastCal);

  return outSoln;
} /* end ObitOTFGetInstCal */

/**
 * Determine instrumental gain from Penn Array-like cal measurements
 * Penn Array-like = slow switching, many samples between state changes.
 * Average On and Off for each detector and differences
 * Mult factor = calJy/(cal_on-cal_off) per detector, may be negative.
 * Data scan averaged and repeated for any subsequent scans without cal On data
 * If doWate==True then determine weights from the reciprocal of the variance
 * in scans which have both cal-on and cal-off measurements; values are reused 
 * between such scans
 * Write Soln entries at beginning and end of each scan.
 * Note: Any scans prior to a scan with calibration data will be flagged.
 * \li "calJy"     OBIT_float (*,1,1) Calibrator value in Jy per detector [def 1.0] .
 *                 Duplicates for all detectors if only one given.
 * \li "doWate"    OBIT_boolean (1,1,1) If true determine Weights from RMS [def False]
 * \li "minSNR"    OBIT_float (1,1,1) Minimum SNR for acceptable cal (ratio cal/RMS)
 *                 Only used if doWate is True [def.  10.0]
 * \param inOTF    Input OTF data, prior calibration/selection applied if requested
 * \param outOTF   OTF with which the output  OTFSoln is to be associated
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created OTFSoln object which is 
 * associated with outOTF.
 */
ObitTableOTFSoln* ObitOTFGetSolnPARGain (ObitOTF *inOTF, ObitOTF *outOTF, ObitErr *err)
{
  ObitTableOTFSoln *outSoln=NULL;
  ObitTableOTFSolnRow *row=NULL;
  ObitOTFDesc *desc=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ObitIOAccess access;
  ofloat *iCal, **accum, *time, blankID;
  ofloat *rec, *calJy=NULL, minsnr, fblank = ObitMagicF();
  ofloat t0, sumTime, lastTarget=-1.0, lastScan=-1.0, lastTime=-1.0;
  ofloat snr, *lastCal=NULL, *lastWeight=NULL;
  gboolean doCalSelect, doWate, fitWate, someData=FALSE, gotCal=FALSE;

  olong iRow, ver, i, j, lrec, ncoef;
  olong  npoly, nDet, nTime, incdatawt;
  ObitIOCode retCode;
  gchar *tname;
  gchar *routine = "ObitOTFGetSolnPARGain";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outSoln;
  g_assert (ObitOTFIsA(inOTF));
  desc = inOTF->myDesc;

  /* Calibration wanted? */ 
  doCalSelect = FALSE;
  ObitInfoListGetTest(inOTF->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

  /* Weights wanted? */ 
  doWate = FALSE;
  ObitInfoListGetTest(inOTF->info, "doWate", &type, dim, &doWate);

  /* If weights are wanted lookup ID of "Blank" scans */
  if (doWate) blankID = FindBlankID(inOTF, err);
  else blankID = -1.0;
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* min SNR */ 
  minsnr = 10.0;
  ObitInfoListGetTest(inOTF->info, "minSNR", &type, dim, &minsnr);

  /* open OTF data to fully instantiate if not already open */
  if ((inOTF->myStatus==OBIT_Inactive) || (inOTF->myStatus==OBIT_Defined)) {
    retCode = ObitOTFOpen (inOTF, access, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
  }

  /* Create output */
  tname   = g_strconcat ("Calibration for: ",inOTF->name, NULL);
  ver     = 0;
  ncoef   = 1;
  npoly   = ncoef;
  outSoln = newObitTableOTFSolnValue(tname, (ObitData*)outOTF, &ver, OBIT_IO_WriteOnly,  
				   inOTF->geom->numberDetect, npoly, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Create work arrays */
  nDet = inOTF->geom->numberDetect;
  lrec = desc->lrec;
  lastCal    = g_malloc0(nDet*sizeof(ofloat)); 
  lastWeight = g_malloc0(nDet*sizeof(ofloat)); 
  calJy      = g_malloc0(nDet*sizeof(ofloat)); 
  iCal       = g_malloc0(MAXSAMPLE*sizeof(ofloat)); 
  time       = g_malloc0(MAXSAMPLE*sizeof(ofloat)); 
  accum      = g_malloc0(nDet*sizeof(ofloat*)); 
  for (i=0; i<nDet; i++) accum[i] = g_malloc0(MAXSAMPLE*sizeof(ofloat));
  for (i=0; i<nDet; i++) lastCal[i]    = fblank;
  for (i=0; i<nDet; i++) lastWeight[i] = fblank;
  t0 = -1.0e20;
  incdatawt = inOTF->myDesc->incdatawt; /* increment in data-wt axis */

  /* cal value in Jy */
  for (j=0; j<nDet; j++) calJy[j] = 1.0;
  ObitInfoListGetTest(inOTF->info, "calJy",  &type, dim, (gpointer*)calJy);
  /* Duplicate if only one */
  if (dim[0]==1) for (j=0; j<nDet; j++) calJy[j] = calJy[0];
 
  /* Open table */
  if ((ObitTableOTFSolnOpen (outSoln, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input OTFSoln table");
    goto cleanup;
  }

  /* Create Row */
  row = newObitTableOTFSolnRow (outSoln);

  /* Attach row to output buffer */
  ObitTableOTFSolnSetRow (outSoln, row, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Initialize solution row */
  row->dAz = 0.0;  /* no pointing corrections */
  row->dEl = 0.0; /* no pointing corrections */
  row->Target = 0;
  for (j=0; j<nDet; j++)  row->mult[j] = 1.0;
  for (j=0; j<nDet; j++)  row->wt[j]   = 1.0;
  for (j=0; j<nDet; j++)  row->cal[j]  = 0.0;
  for (j=0; j<nDet; j++)  row->add[j]  = 0.0;
  for (i=0; i<npoly; i++) row->poly[i] = 0.0;

  /* loop over data file */
  retCode = OBIT_IO_OK;
  sumTime = 0.0;
  nTime = 0;
  while (retCode == OBIT_IO_OK) {
    
    /* read buffer */
    if (doCalSelect) retCode = ObitOTFReadSelect (inOTF, NULL, err);
    else retCode = ObitOTFRead (inOTF, NULL, err);
    if (err->error) goto cleanup;
    if (retCode==OBIT_IO_EOF) break; /* done? */
    
    /* Record pointer */
    rec = inOTF->buffer;
    
    /* First time */
    if (t0<-1.0e10) {
      t0 = rec[desc->iloct];
      lastTime   = rec[desc->iloct];
      lastTarget = rec[desc->iloctar];
      lastScan   = rec[desc->ilocscan];
    }
    
    /* Loop over buffer */
    for (i=0; i<desc->numRecBuff; i++) {
      
      /* Accumulation (scan) finished? If so, compute solution and write.*/
      if ((rec[desc->iloctar]  != lastTarget) || (nTime>=MAXSAMPLE) ||
	  (rec[desc->ilocscan] != lastScan)) {
	
	/* Not first time - assume first descriptive parameter never blanked */
	if (nTime>0) {
	  /* Set descriptive info on Row */
	  row->Time  = sumTime/nTime;  /* center time */
	  row->TimeI = 2.0 * (row->Time - t0);

	  /* Flatten curves, estimate cal, Weights */
	  /* Is this a "Blank Scan? */
	  fitWate = (fabs(lastTarget-blankID)<0.1) || (blankID<0.0);
	  doPARCalWeight(nDet, nTime, iCal, accum, time, lastCal, lastWeight, fitWate);
	  
	  /* Propagate last good cals if needed */
	  for (j=0; j<nDet; j++ ) {
	    row->cal[j]  = lastCal[j];
	    if (lastCal[j]!=fblank) row->mult[j] = calJy[j] / lastCal[j];
	    else row->mult[j] = 1.0;
	  }

	  /* Calibrating weights? Propagate last weight */
	  if (doWate) {
	    for (j=0; j<nDet; j++ ) {
	      /* Acceptable SNR? */
	      if ((lastCal[j]!=fblank) && (lastWeight[j]!=fblank)) {
		snr = fabs(lastCal[j]) * sqrt(lastWeight[j]);
		if (snr < minsnr) {
		  /* No - off with their heads! */
		  row->wt[j]  = 0.0;
		  row->add[j] = fblank;
		}
	      }  /* End check SNR */
	      if (lastWeight[j]!=fblank) {
		row->wt[j]  = lastWeight[j];
		row->add[j] = 0.0;
	      }
	    } /* end loop over detectors */
	  } /* end calibrating weights */
	  
	  /* Write Soln table - bracket interval */
	  iRow = -1;
	  row->Target = (oint)(lastTarget+0.5);
	  row->Time   = t0;
	  if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	       != OBIT_IO_OK) || (err->error>0)) { 
	    Obit_log_error(err, OBIT_Error, "ERROR writing OTFSoln Table file");
	    goto cleanup;
	  }
	  row->Time   = lastTime;
	  if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	       != OBIT_IO_OK) || (err->error>0)) { 
	    Obit_log_error(err, OBIT_Error, "ERROR writing OTFSoln Table file");
	    goto cleanup;
	  }
	  /* initialize accumulators */
	  t0 = rec[desc->iloct];
	  lastTime   = rec[desc->iloct];
	  lastTarget = rec[desc->iloctar];
	  lastScan   = rec[desc->ilocscan];
	  sumTime = 0.0;
	  nTime = 0;
	} /* end of do solution if there is data */
	
      } /* end do solution */
      
      /* accumulate */
      lastTime = rec[desc->iloct];
      sumTime += rec[desc->iloct];
      iCal[nTime] = fabs(rec[inOTF->myDesc->iloccal]);
      time[nTime] = rec[inOTF->myDesc->iloct];
      for (j=0; j<nDet; j++ ) {
	accum[j][nTime] = rec[desc->ilocdata+j*incdatawt];
      }
      nTime++;
      rec += desc->lrec; /* Data record pointer */
      someData = TRUE;
      
    } /* end loop over buffer load */
  } /* end loop reading data file */
  
  /* Finish up any data in accumulator */
  if (nTime>0) {
    /* Set descriptive info on Row */
    row->Time  = sumTime/nTime; /* time */
    row->TimeI = 2.0 * (row->Time - t0);
    
    /* Flatten curves, estimate cal, Weights */
    fitWate = (fabs(lastTarget-blankID)<0.1) || (blankID<0.0);
    doPARCalWeight(nDet, nTime, iCal, accum, time, lastCal, lastWeight, fitWate);
	  
    /* Propagate last good cals if needed */
    for (j=0; j<nDet; j++ ) {
      row->cal[j]  = lastCal[j];
      if (lastCal[j]!=fblank) row->mult[j] = calJy[j] / lastCal[j];
      else row->mult[j] = 1.0;
    }

    /* Calibrating weights? Propagate last weight */
    if (doWate) {
      for (j=0; j<nDet; j++ ) {
	/* Acceptable SNR? */
	if ((lastCal[j]!=fblank) && (lastWeight[j]!=fblank)) {
	  snr = fabs(lastCal[j]) * sqrt(lastWeight[j]);
	  if (snr < minsnr) {
	    /* No - off with their heads! */
	    row->wt[j]  = 0.0;
	    row->add[j] = fblank;
	  }
	}  /* End check SNR */
	if (lastWeight[j]!=fblank) row->wt[j] = lastWeight[j];
      } /* end loop over detectors */
    } /* end calibrating weights */
	  
    /* Write Soln table - bracket interval */
    iRow = -1;
    row->Time = t0;
    if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR writing OTFSoln Table file");
      goto cleanup;
    }
    row->Time   = lastTime;
    if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR writing OTFSoln Table file");
      goto cleanup;
    }
  } /* End final solution */

  /* Make sure this worked -final lastCal must have some good data */
  gotCal = FALSE;
  for (j=0; j<nDet; j++ ) {
    if ((lastCal[j]!=fblank) && (lastWeight[j]!=fblank)) gotCal = TRUE;
  }
  if (!gotCal) {
    Obit_log_error(err, OBIT_Error, "%s: ERROR NO calibration data found", routine);
  }

  /* Cleanup */
 cleanup:
  /* Close cal table */
  if ((ObitTableOTFSolnClose (outSoln, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input OTFSoln Table file");
  }
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);

  /* Give warning if no data selected */
  if (!someData) Obit_log_error(err, OBIT_InfoWarn, 
				"%s: Warning: NO data selected", routine);

  /* Cleanup */
  if (calJy)      g_free(calJy);
  if (lastCal)    g_free(lastCal);
  if (lastWeight) g_free(lastWeight);
  for (i=0; i<nDet; i++) if (accum[i]) g_free(accum[i]);
  if (accum)   g_free(accum);
  if (iCal)    g_free(iCal);
  if (time)    g_free(time);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  return outSoln;
} /* end ObitOTFGetSolnPARGain */

/**
 * Write pointing corrections into a table
 * Given a set of pointing correction, write to OTFSoln table
 * Calibration parameters are on the inOTF info member.
 * \li "POffset"   OBIT_float (3,?,1) Table of pointing offsets
 *                 Triplets (time(day), Azoff(asec), Decoff(asec))
 * \param inOTF    Input OTF data. 
 * \param outOTF   OTF with which the output OTFSoln is to be associated
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created OTFSoln object which is 
 * associated with inOTF.
 */
ObitTableOTFSoln* ObitOTFGetSolnPointTab (ObitOTF *inOTF, ObitOTF *outOTF, ObitErr *err)
{
  ObitTableOTFSoln *outSoln=NULL;
  ObitTableOTFSolnRow *row=NULL;
  gint32 pntdim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat *pointTab;
  ObitInfoType type;
  olong iRow, ver, j, mpoly, ndetect;
  gchar *tname;
  gchar *routine = "ObitOTFGetSolnPointTab";
 
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outSoln;
  g_assert (ObitOTFIsA(inOTF));

  /* Create output */
  tname = g_strconcat ("Calibration for: ",inOTF->name, NULL);
  ver = 0;
  mpoly = 1;
  outSoln = newObitTableOTFSolnValue(tname, (ObitData*)outOTF, &ver, OBIT_IO_WriteOnly,  
				     inOTF->geom->numberDetect, mpoly, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Get offset table  */
  pointTab = NULL;
  ObitInfoListGetP(inOTF->info, "POffset",   &type, pntdim, (gpointer*)&pointTab);
  /* error checks */
  if (pointTab==NULL) {
    Obit_log_error(err, OBIT_Error, "%s ERROR, NO POffset array given", routine);
    goto cleanup;
  }
  if ((pntdim[0]!=3) || (type!=OBIT_float)) {
    Obit_log_error(err, OBIT_Error, "%s Invalid POffset array given", routine);
    goto cleanup;
  }

  /* Open output table */
  if ((ObitTableOTFSolnOpen (outSoln, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "%s ERROR opening input OTFSoln table", routine);
    goto cleanup;
  }

  /* Create Row */
  row = newObitTableOTFSolnRow (outSoln);

  /* Attach row to output buffer */
  ObitTableOTFSolnSetRow (outSoln, row, err);
  if (err->error) goto cleanup;

  /* Initialize solution row */
  ndetect = inOTF->geom->numberDetect;  /* number of detectors */
  row->dAz = 0.0;  /* no pointing corrections */
  row->dEl = 0.0; /* no pointing corrections */
  row->Target = 0;
  row->TimeI = 0.0;
  for (j=0; j<ndetect; j++) row->mult[j] = 1.0;
  for (j=0; j<ndetect; j++) row->wt[j]   = 1.0;
  for (j=0; j<ndetect; j++) row->cal[j]  = 0.0;
  for (j=0; j<ndetect; j++) row->add[j]  = 0.0;

  /* Loop over table */
  for (iRow=1; iRow<=pntdim[1]; iRow++) {
    row->Time = pointTab[iRow*3-3];
    row->dAz  = pointTab[iRow*3-2]/3600.0;  /* To deg */
    row->dEl  = pointTab[iRow*3-1]/3600.0;
    if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "%s ERROR writing OTFSoln Table file", routine);
      goto cleanup;
    }
  } /* end loop writing table */

  /* Close output cal table */
  if ((ObitTableOTFSolnClose (outSoln, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s ERROR closing output OTFSoln Table file", routine);
    goto cleanup;
  }
  
  /* Deallocate arrays */
  cleanup:
  row = ObitTableOTFSolnUnref(row);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  return outSoln;
} /* end ObitOTFGetSolnPointTab */

/**
 * Create dummy OTFCal table table (applying will not modify data)
 * \li "solInt"    OBIT_float (1,1,1) Solution interval in days [def 1 sec].
 * \param inOTF    Input OTF data. 
 * \param outOTF   OTF with which the output  OTFCal is to be associated
 * \param ver      OTFCal table version
 * \param ncoef    Number of coefficients in table
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created OTFCal object which is 
 * associated with inOTF.
 */
ObitTableOTFCal* ObitOTFGetDummyCal (ObitOTF *inOTF, ObitOTF *outOTF, olong ver, 
				     oint ncoef, ObitErr *err)
{
#define MAXSAMPLE 1000   /* Maximum number of samples in an integration */
  ObitTableOTFCal *outCal=NULL;
  ObitTableOTFCalRow *row=NULL;
  ObitOTFDesc *desc=NULL;
  ObitIOAccess access;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat *rec, solInt, t0, sumTime, lastTime=-1.0, lastTarget=-1.0, lastScan=-1.0;
  olong iRow, i, lrec;
  olong  nDet, nTime;
  gboolean doCalSelect, someData=FALSE;
  ObitIOCode retCode;
  gchar *tname;
  gchar *routine = "ObitOTFGetDummyCal";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outCal;
  g_assert (ObitOTFIsA(inOTF));
  desc = inOTF->myDesc;

  /* Calibration/selection wanted? */ 
  doCalSelect = FALSE;
  ObitInfoListGetTest(inOTF->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

  /* open OTF data to fully instantiate if not already open */
  if ((inOTF->myStatus==OBIT_Inactive) || (inOTF->myStatus==OBIT_Defined)) {
    retCode = ObitOTFOpen (inOTF, access, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);
  }
  lrec = inOTF->myDesc->lrec;
  nDet = inOTF->geom->numberDetect;
  t0 = -1.0e20;

  /* Create output */
  tname = g_strconcat ("Calibration for: ",inOTF->name, NULL);
  outCal = newObitTableOTFCalValue(tname, (ObitData*)outOTF, &ver, OBIT_IO_WriteOnly,  
				   inOTF->geom->numberDetect, ncoef, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);

  /* Get parameters for calibration */
  /* "Solution interval" default 1 sec */
  solInt = 1.0 / 86400.0;
  ObitInfoListGetTest(inOTF->info, "solInt", &type, dim, (gpointer*)&solInt);

  /* Open table */
  if ((ObitTableOTFCalOpen (outCal, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening input OTFCal table", routine);
    return outCal;
  }

  /* Create Row */
  row = newObitTableOTFCalRow (outCal);

  /* Attach row to output buffer */
  ObitTableOTFCalSetRow (outCal, row, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);

  /* Initialize */
  row->dAz = 0.0;  /* no pointing corrections */
  row->dEl = 0.0; /* no pointing corrections */
  for (i=0; i<nDet; i++) {
    row->mult[i] = 1.0;
    row->wt[i]   = 1.0;
    row->add[i]  = 0.0;
    row->cal[i]  = 0.0;
  }

  /* loop calibrating data */
  retCode = OBIT_IO_OK;
  sumTime = 0.0;
  nTime = 0;
  while (retCode == OBIT_IO_OK) {
    
    /* read buffer */
    if (doCalSelect) retCode = ObitOTFReadSelect (inOTF, NULL, err);
    else retCode = ObitOTFRead (inOTF, NULL, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);
    if (retCode==OBIT_IO_EOF) break; /* done? */
    
    /* Record pointer */
    rec = inOTF->buffer;
    
    /* First time */
    if (t0<-1.0e10) {
      t0         = rec[inOTF->myDesc->iloct];
      lastTarget = rec[inOTF->myDesc->iloctar];
      lastScan   = rec[inOTF->myDesc->ilocscan];
      lastTime   = rec[inOTF->myDesc->iloct];
      /* Write entry at very beginning */
      row->Time  = rec[inOTF->myDesc->iloct]; 
      row->TimeI = 0.0;
      row->Target = (oint)(rec[inOTF->myDesc->iloctar]+0.5);
      iRow = -1;
      if ((ObitTableOTFCalWriteRow (outCal, iRow, row, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "%s: ERROR writing OTFCal Table file", routine);
	return outCal;
      }
    }
    
    /* Loop over buffer */
    for (i=0; i<inOTF->myDesc->numRecBuff; i++) {

      /* Accumulation or scan finished? If so, write "calibration".*/
      if ((rec[inOTF->myDesc->iloct] > (t0+solInt)) || (nTime>=MAXSAMPLE) ||
	  (rec[inOTF->myDesc->iloctar] != lastTarget) ||  
	  (rec[inOTF->myDesc->ilocscan] != lastScan)) {
	
	/* Not first time - assume first descriptive parameter never blanked */
	if (nTime>0) {
	  /* if new scan write end of last scan and this time */
	  if ((rec[inOTF->myDesc->iloctar] != lastTarget) ||  
	      (rec[inOTF->myDesc->ilocscan] != lastScan)) {

	    /* values for end of previous scan */
	    row->Time  = lastTime; 
	    row->TimeI = 0.0;
	    row->Target = (oint)(lastTarget+0.5);
	    iRow = -1;
	    if ((ObitTableOTFCalWriteRow (outCal, iRow, row, err)
		 != OBIT_IO_OK) || (err->error>0)) { 
	      Obit_log_error(err, OBIT_Error, "%s: ERROR writing OTFCal Table file", routine);
	      return outCal;
	    }
	    /* Values for start of next scan */
	    row->Time  = rec[inOTF->myDesc->iloct]; 
	    row->TimeI = 0.0;
	    row->Target = (oint)(rec[inOTF->myDesc->iloctar]+0.5);
	  } else {  /* in middle of scan - use average time */
	    /* Set descriptive info on Row */
	    row->Time  = sumTime/nTime;  /* time */
	    row->TimeI = 2.0 * (row->Time - t0);
	    row->Target = (oint)(rec[inOTF->myDesc->iloctar]+0.5);
	  }
      
	  /* Write Cal table */
	  iRow = -1;
	  if ((ObitTableOTFCalWriteRow (outCal, iRow, row, err)
	       != OBIT_IO_OK) || (err->error>0)) { 
	    Obit_log_error(err, OBIT_Error, "%s: ERROR writing OTFCal Table file", routine);
	    return outCal;
	  }
	  /* initialize accumulators */
	  t0         = rec[inOTF->myDesc->iloct];
	  lastTarget = rec[inOTF->myDesc->iloctar];
	  lastScan   = rec[inOTF->myDesc->ilocscan];
	  sumTime = 0.0;
	  nTime = 0;
	} /* end of do solution if there is data */
	
      } /* end do solution */
      
      /* accumulate */
      sumTime += rec[inOTF->myDesc->iloct];
      lastTime = rec[inOTF->myDesc->iloct];
      nTime++; /* how many data points */
      rec += inOTF->myDesc->lrec; /* Data record pointer */
      someData = TRUE;
      
    } /* end loop over buffer load */
  } /* end loop reading/gridding data */
  
  /* Finish up any data in accumulator */
  if (nTime>0) {
    /* Set descriptive info on Row */
    row->Time  = lastTime; 
    row->TimeI = 0.0;
    row->Target = (oint)(lastTarget+0.5);
    /* Write Cal table */
    iRow = -1;
    if ((ObitTableOTFCalWriteRow (outCal, iRow, row, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "%s: ERROR writing OTFCal Table file", routine);
      return outCal;
    }
  } /* End final cal */

  /* Close cal table */
  if ((ObitTableOTFCalClose (outCal, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing input OTFCal Table file", routine);
    return outCal;
  }
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);

  /* Give warning if no data selected */
  if (!someData) Obit_log_error(err, OBIT_InfoWarn, 
				"%s: Warning: NO data selected", routine);

  return outCal;
} /* end ObitOTFGetDummyCal */

/**
 * Determine flagging from the statistics of time stream data.
 * If a model is supplied, this is determined from the residuals to the model.
 * The residual data will be derived from inOTF minus the model.
 * Detector/intervals in excess of the larger of maxRMS or maxRatio times the 
 * equivalent model RMS are flagged in OTFFlag table FGVer on outOTF.
 * An evaluation is made independently in each flagInt and detector.
 * Control parameters are on the inOTF info member.
 * \li "flagInt"   OBIT_float (1,1,1) Flaffing interval in days [def 10 sec].
 *                 This should not exceed 1000 samples.  Intervals will be truncated
 *                 at this limit.
 * \li "maxRMS"    OBIT_float (1,1,1) Maximum allowable  RMS in Jy[ def 1.0].
 * \li "maxRatio"  OBIT_float (*,1,1) Max. allowable ratio to equivalent model RMS [2]
 * \li "minEl"     OBIT_float (1,1,1) Minimum elevation allowed (deg)
 * \param inOTF    Input OTF data. 
 * \param model    Input CLEAN model, if not provided (NULL) then the model RMS 
 *                 is assumed 0.0
 *                 The pixels values in the image array should be the estimated
 *                 antenna response to the source in that direction.
 * \param outOTF   OTF with which the output OTFFlag is to be associated
 * \param FGVer    Flag version for output, if 0 create new
 * \param err      Error stack, returns if not empty.
 */
void ObitOTFGetSolnFlag (ObitOTF *inOTF, ObitImage *model, 
			 ObitOTF *outOTF, olong FGVer, ObitErr *err)
{
#define MAXFLAGSAMPLE 1000   /* Maximum number of samples in an integration */
  ObitTableOTFFlag *outFlag=NULL, *inFlag=NULL;
  ObitTableOTFFlagRow *row=NULL;
  ObitOTF *modelOTF=NULL, *residOTF=NULL;
  ObitIOAccess access;
  ObitIOSize IOBy;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat lastScan=-1.0, lastTarget=-1.0;
  ofloat *residRec=NULL, *modelRec=NULL, flagInt, minEl, maxRMS, maxRatio, el, t0, tEnd;
  ofloat *time, *data=NULL, *moddata=NULL, residRMS, modelRMS;
  ofloat fblank = ObitMagicF();
  olong iRow, oRow,ver, i, j, lrec, nsample, flagCnt, totalCnt;
  olong  ndetect, incdatawt, NPIO;
  ObitIOCode retCode;
  gboolean flag, someData, doCalSelect, doFlag;
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar *tname;
  gchar *routine = "ObitOTFGetSolnFlag";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return ;
  g_assert (ObitOTFIsA(inOTF));
  g_assert (ObitOTFIsA(outOTF));

  /* Create output - write new OTFFlag table and then copy to FGVer at the end
     This avoids the potential for trouble of reading (applying) and writing the 
     same table at the same time */
  tname = g_strconcat ("Flagging for: ",outOTF->name, NULL);
  ver = 0;
  outFlag = newObitTableOTFFlagValue(tname, (ObitData*)outOTF, &ver, OBIT_IO_WriteOnly,  
				     err);
  g_free (tname);
  if (err->error) Obit_traceback_msg (err, routine, outOTF->name);

  /* If model given generate model and residual data files */
  NPIO = 1000;
  dim[0] = 1;
  ObitInfoListAlwaysPut (inOTF->info, "nRecPIO", OBIT_long, dim,  &NPIO);
  if (model) {
    /* Open model */
    IOBy = OBIT_IO_byPlane;
    dim[0] = 1;
    ObitInfoListAlwaysPut (model->info, "IOBy", OBIT_long, dim, &IOBy);
    dim[0] = 7;
    ObitInfoListAlwaysPut (model->info, "BLC", OBIT_long, dim, blc); 
    ObitInfoListAlwaysPut (model->info, "TRC", OBIT_long, dim, trc);
    ObitImageOpen (model, OBIT_IO_ReadOnly, err);
    /* Read input plane */
    ObitImageRead (model, NULL , err);

    modelOTF = newObitOTFScratch(inOTF, err);  /* Scratch file for model */
    dim[0] = 1;
    ObitInfoListPut (modelOTF->info, "nRecPIO", OBIT_long, dim,  &NPIO, err);
    /* Form model data */
    ObitOTFUtilModelImage (inOTF, modelOTF, model->image, model->myDesc, err);
    /* open model OTF data to fully instantiate if not already open */
    ObitOTFOpen (modelOTF, OBIT_IO_ReadOnly, err);
    if (err->error) Obit_traceback_msg (err, routine, modelOTF->name);

    residOTF = newObitOTFScratch(inOTF, err); /* Scratch file for residual */
    /* Form residual data set */
    ObitOTFUtilSubImage (inOTF, residOTF, model->image, model->myDesc, err);

    /* Close model */
    ObitImageClose(model, err);
    if (err->error) Obit_traceback_msg (err, routine, residOTF->name);
    /* free image memory if not memory resident */
    if (model->mySel->FileType!=OBIT_IO_MEM) 
      model->image = ObitFArrayUnref(model->image);
  } else { /* No model - use data as residual */
    residOTF = ObitOTFRef(inOTF);
  }

  /* Calibration wanted for residual? (only really if no model) */ 
  doCalSelect = FALSE;
  ObitInfoListGetTest(residOTF->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

  /* open resid OTF data */
  dim[0] = 1;
  ObitInfoListPut (residOTF->info, "nRecPIO", OBIT_long, dim,  (gpointer)&NPIO, err);
  ObitOTFOpen (residOTF, access, err);
  if (err->error) Obit_traceback_msg (err, routine, residOTF->name);

  /* Create work arrays */
  ndetect  = inOTF->geom->numberDetect;  /* number of detectors */
  lrec     = inOTF->myDesc->lrec;
  time     = g_malloc0(MAXFLAGSAMPLE*sizeof(ofloat));
  if (model) moddata  = g_malloc0(ndetect*MAXFLAGSAMPLE*sizeof(ofloat));
  data     = g_malloc0(ndetect*MAXFLAGSAMPLE*sizeof(ofloat));
  nsample  = 0;
  t0       = -1.0e20;
  tEnd     = -1.0e20;
  lastScan = -1000;
  incdatawt= inOTF->myDesc->incdatawt; /* increment on data-wt axis */
  someData = FALSE;

  /* Get parameters for flagging */
  /* interval default 10 sec */
  flagInt = 10.0 / 86400.0;
  ObitInfoListGetTest(inOTF->info, "flagInt", &type, dim, (gpointer*)&flagInt);
  Obit_log_error(err, OBIT_InfoErr, 
		 "%s:  Flagging interval %6.2f second", routine, flagInt*86400.0);
  
  /* maximum RMS in interval */
  maxRMS = 1.0;
  ObitInfoListGetTest(inOTF->info, "maxRMS", &type, dim, (gpointer*)&maxRMS);

  /* maximum ratio to model RMS in interval */
  maxRatio = 2.0;
  ObitInfoListGetTest(inOTF->info, "maxRatio", &type, dim, (gpointer*)&maxRatio);

  /* minimum allowed elevation for interval */
  minEl = 1.0;
  ObitInfoListGetTest(inOTF->info, "minEl",  &type, dim, (gpointer*)&minEl);

  /* Open output table */
  if ((ObitTableOTFFlagOpen (outFlag, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input OTFFlag table");
    goto cleanup;
  }

  /* Create Row */
  row = newObitTableOTFFlagRow (outFlag);

  /* Attach row to output buffer */
  ObitTableOTFFlagSetRow (outFlag, row, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Initialize solution row */
  row->TargetID = 0;
  row->Feed = (oint)0;
  row->TimeRange[0] = 0.0; row->TimeRange[1] = 0.0; 
  row->chans[0] = 1; row->chans[1] = 0;
  row->pFlags[0] = 0x1 | 0x2 | 0x4 | 0x8;  /* Flag all Stokes */
  strcpy (row->reason, "Excessive RMS");
  flag = FALSE;
  residRec = residOTF->buffer;
  if (model) modelRec = modelOTF->buffer;
  flagCnt = 0;
  totalCnt = 0;

  /* loop over data */
  retCode = OBIT_IO_OK;
  while (retCode == OBIT_IO_OK) {

    /* read buffer */
    if (doCalSelect) retCode = ObitOTFReadSelect (residOTF, NULL, err);
    else retCode = ObitOTFRead (residOTF, NULL, err);
    if (err->error) goto cleanup;
    if (retCode==OBIT_IO_EOF) break; /* done? */
    /* Read model if given */
    if (model) {
      retCode = ObitOTFRead (modelOTF, NULL, err);
      if (err->error) goto cleanup;
      if (retCode==OBIT_IO_EOF) break; /* done? */
    }

    /* Record pointers */
    residRec = residOTF->buffer;
    if (model) {
      modelRec = modelOTF->buffer;
      /* Also consistency checks - check number of records in buffer */
      if (modelOTF->myDesc->numRecBuff!=residOTF->myDesc->numRecBuff) {
	Obit_log_error(err, OBIT_Error, 
		       "%s: model and residual have different no. records", 
		       routine);
	goto cleanup;
      }
      /* Check timetags of first records */
      if (residRec[residOTF->myDesc->iloct]!=modelRec[modelOTF->myDesc->iloct]) {
	Obit_log_error(err, OBIT_Error, 
		       "%s: model and residual have different timetags", 
		       routine);
	goto cleanup;
      }
    } /* End consistency checks */
  
    /* First time */
    if (t0<-1.0e10) {
      t0 = residRec[residOTF->myDesc->iloct]; 
      tEnd = t0+1.0E-20;
      lastScan   = residRec[residOTF->myDesc->ilocscan]; /* Which scan number */
      lastTarget = residRec[residOTF->myDesc->iloctar];  /* Which target number */
    }

    /* Loop over buffer */
    for (i=0; i<residOTF->myDesc->numRecBuff; i++) {

      /* Flagging interval finished? If so, compute solution and write.
	 also if work arrays full */
      if ((residRec[residOTF->myDesc->iloct] > (t0+flagInt)) || 
	  /* or end of scan */
	  ( residRec[residOTF->myDesc->ilocscan] != lastScan) ||
	  /* or exceed limit on number of samples */
	  (nsample>=MAXFLAGSAMPLE)) {
	
	/* Any data in buffers? */
	if (nsample>0) {
	/* Set descriptive info on Row */
	  row->TargetID = (oint)(lastTarget+0.5);
	  row->TimeRange[0] = t0;  /* time range */
	  row->TimeRange[1] = tEnd;
	  
	  /* Check each detector */
	  doFlag = FALSE;
	  /* Is data below the elevation limit? */
	  el = ObitOTFArrayGeomElev (residOTF->geom, residRec[residOTF->myDesc->iloct], 
				     residRec[residOTF->myDesc->ilocra], 
				     residRec[residOTF->myDesc->ilocdec]);
	  for (j=0; j<ndetect; j++) {
	    /* RMSes of residual, model */
	    residRMS = RMSValue (&data[j*MAXFLAGSAMPLE], nsample);
	    if (residRMS!=fblank) {
	      totalCnt++;
	      if (model) {
		/* have model */
		modelRMS = RMSValue (&moddata[j*MAXFLAGSAMPLE], nsample);
		if (modelRMS!=fblank)  doFlag = residRMS > MAX(maxRatio*modelRMS, maxRMS);
	      } else {
		/* Only residual */
		doFlag = residRMS > maxRMS;
	      }

	      doFlag = doFlag || (el < minEl);
	    }
	    
	    /* Write entry if needed */
	    if (doFlag) {
	      row->Feed = (oint)(j+1);
	      iRow = -1;
	      ObitTableOTFFlagWriteRow (outFlag, iRow, row, err);
	      if (err->error) goto cleanup;
	      flagCnt++;
	    }
	  } /* End loop checking detectors */
	} /* end if some valid data */
	

	/* initialize accumulators */
	t0 = residRec[residOTF->myDesc->iloct];
	tEnd = t0+1.0E-20;
	nsample = 0;
	
      } /* end do solution */
      
      /* accumulate */
      if (!flag) {
	time[nsample] = residRec[residOTF->myDesc->iloct];
	for (j=0; j<ndetect; j++ ) 
	  data[nsample+j*MAXFLAGSAMPLE] = residRec[residOTF->myDesc->ilocdata+j*incdatawt];
	if (model) {
	  for (j=0; j<ndetect; j++ ) 
	    moddata[nsample+j*MAXFLAGSAMPLE] = modelRec[modelOTF->myDesc->ilocdata+j*incdatawt];
	}
	nsample++;
	someData = TRUE;  /* Any valid data? */
      }
      /* last time in accumulation */
      tEnd       = residRec[residOTF->myDesc->iloct]+1.0E-20;
      lastScan   = residRec[residOTF->myDesc->ilocscan]; /* Which scan number */
      lastTarget = residRec[residOTF->myDesc->iloctar];  /* Which target number */
      residRec += residOTF->myDesc->lrec;            /* Data record pointers */
      if (model) modelRec += modelOTF->myDesc->lrec; 
      
    } /* end loop over buffer load */
  } /* end loop reading data */
  
  /* Finish up any data in arrays */
  /* Any data in buffers? */
  if (nsample>0) {
    /* Set descriptive info on Row */
    row->TargetID = (oint)(lastTarget+0.5);
    row->TimeRange[0] = t0;  /* time range */
    row->TimeRange[1] = tEnd;
    
    /* Check each detector */
    doFlag = FALSE;
    for (j=0; j<ndetect; j++) {
      /* RMSes of residual, model */
      residRMS = RMSValue (&data[j*MAXFLAGSAMPLE], nsample);
      if (residRMS!=fblank) {
	totalCnt++;
	if (model) {
	  /* have model */
	  modelRMS = RMSValue (&moddata[j*MAXFLAGSAMPLE], nsample);
	  if (modelRMS!=fblank)  doFlag = residRMS > MAX(maxRatio*modelRMS, maxRMS);
	} else {
	  /* Only residual */
	  doFlag = residRMS > maxRMS;
	}
      }

      /* Is data below the elevation limit? */
      el = ObitOTFArrayGeomElev (residOTF->geom, residRec[residOTF->myDesc->iloct], 
				 residRec[residOTF->myDesc->ilocra], 
				 residRec[residOTF->myDesc->ilocdec]);
      doFlag = doFlag || (el < minEl);
      
      /* Write entry if needed */
      if (doFlag) {
	row->Feed = (oint)(j+1);
	iRow = -1;
	ObitTableOTFFlagWriteRow (outFlag, iRow, row, err);
	if (err->error) goto cleanup;
	flagCnt++;
      }
    } /* End loop checking detectors */
  } /* end if some valid data */
  
  /* tell how much flagged */
  Obit_log_error(err, OBIT_InfoErr, 
		 "%s:  Flagged  %d/ %d detector/intervals", routine, flagCnt, totalCnt);
  
  /* Give warning if no data selected */
  if (!someData) Obit_log_error(err, OBIT_InfoWarn, 
				"%s: Warning: NO data selected", routine);

  /* if ver is not FGVer, then copy entries and delete the temporary table */
  if (FGVer<=0) FGVer = ver;
  if (FGVer!=ver) {
    /* Clean up old */
    row      = ObitTableOTFFlagUnref(row);
    outFlag  = ObitTableOTFFlagUnref(outFlag);

    /* Input table */
    tname = g_strconcat ("input for: ",outOTF->name, NULL);
    inFlag = newObitTableOTFFlagValue(tname, (ObitData*)outOTF, &ver, OBIT_IO_ReadOnly,  
				       err);
    g_free (tname);
    if (err->error) Obit_traceback_msg (err, routine, outOTF->name);

    /* Output table */
    tname = g_strconcat ("output for: ",outOTF->name, NULL);
    outFlag = newObitTableOTFFlagValue(tname, (ObitData*)outOTF, &FGVer, OBIT_IO_ReadWrite,  
				       err);
    g_free (tname);
    if (err->error) Obit_traceback_msg (err, routine, outOTF->name);

    /* append input to output */
    /* Open input tables */
    ObitTableOTFFlagOpen (inFlag, OBIT_IO_ReadOnly, err);
    ObitTableOTFFlagOpen (outFlag, OBIT_IO_ReadWrite, err);
    if (err->error) goto cleanup;
    
    /* If there are entries in the output table, mark it unsorted */
    if (outFlag->myDesc->nrow>0) 
      {outFlag->myDesc->sort[0]=0; outFlag->myDesc->sort[1]=0;}
  
    /* Create Row */
    row = newObitTableOTFFlagRow (outFlag);
    
    /* Attach row to output buffer */
    ObitTableOTFFlagSetRow (outFlag, row, err);
    if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

    for (iRow=1; iRow<=inFlag->myDesc->nrow; iRow++) {
      ObitTableOTFFlagReadRow (inFlag, iRow, row, err);
      oRow = -1;
      ObitTableOTFFlagWriteRow (outFlag, oRow, row, err);
      if (err->error) goto cleanup;
    } /* end loop copying */
    
    ObitTableOTFFlagClose (inFlag, err);
    ObitTableOTFFlagClose (outFlag, err);
    /* Zap temporary table */
    ObitOTFZapTable(outOTF, "OTFFlag", ver, err);
    if (err->error) goto cleanup;
  } /* end copy table */

  /* Cleanup */
 cleanup:
  /* Close data */
  retCode = ObitOTFClose (residOTF, err);

  /* Close flag table */
  if ((ObitTableOTFFlagClose (outFlag, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s ERROR closing output OTFFlag Table file", routine);
  }
  
  outFlag  = ObitTableOTFFlagUnref(outFlag);
  inFlag   = ObitTableOTFFlagUnref(inFlag);
  modelOTF = ObitOTFUnref(modelOTF);
  residOTF = ObitOTFUnref(residOTF);
  row      = ObitTableOTFFlagUnref(row);
  if (time)    g_free(time);
  if (data)    g_free(data);
  if (moddata) g_free(moddata);
  if (err->error) Obit_traceback_msg (err, routine, residOTF->name);

} /* end ObitOTFGetSolnFlag */

/*----------------------Private functions---------------------------*/

/**
 * Fit an atmospheric model across a grid of detectors.
 * If a model fitter is provided, it is used, else a median of the values.
 * \param geom   Array geometry structure
 * \param desc   OTF descriptor
 * \param rec    OTF Data record
 * \param fitter Atmospheric model fitter, NULL if none
 * \param row    Solution table row
 */
static void 
ObitOTFGetSolnSolve (ObitOTFArrayGeom *geom, ObitOTFDesc *desc, 
		    ObitPennArrayAtmFit *fitter, ofloat *rec, ObitTableOTFSolnRow *row)
{
  olong i, incdatawt;
  ofloat value, corr, fblank = ObitMagicF();

  /* Initialize */
  row->dAz  = 0.0;  /* no pointing corrections */
  row->dEl = 0.0;  /* no pointing corrections */
  incdatawt = desc->incdatawt; /* increment in data-wt axis */

  /* Model fitting or median? */
  if (fitter!=NULL) { /* model fitting */
    ObitPennArrayAtmFitFit (fitter, &rec[desc->ilocdata], incdatawt, row->poly);
    corr = 0.0;
  } else { /* Get median */
    value = GetSolnMedianAvg (&rec[desc->ilocdata], incdatawt, geom->numberDetect);
    if (value!=fblank) corr = -value;
    else corr = fblank;
  }

   /* Set values on calibration row */
   for (i=0; i<geom->numberDetect; i++) {
     row->add[i]  = corr;
     row->mult[i] = 1.0;
     row->wt[i]   = 1.0;
   }
   
} /* end ObitOTFGetSolnSolve */

/**
 * Use cal measurements in a data set to estimate the gain calibration factor.
 * The value in cal is assumed 0.0 if the cal is off and the value of the cal in Jy if > 0.
 * If a model fitter is provided, it is used, else a median of the values.
 * The average noise cal values are written on row->cal[detect].
 * \param nsample  The number of samples
 * \param time     Array of time tags, used to interpolate cal off values to time of cal on.
 * \param cal      >0 if cal on else cal is off per integration.
 * \param calJy    Cal value in Jy.
 * \param data     data values, per integration
 * \param minrms   minimum allowed RMS in solution as a fraction of the gain.
 * \param lastgood Last good value of gain.
 * \param detect   The detector number (used to determine which value is set in row).
 * \param row      Solution table row
 */
static void 
ObitOTFGetSolnGainSolve (olong nsample, ofloat *time, ofloat *cal, ofloat calJy, ofloat *data, 
			 ofloat minrms, ofloat lastgood, olong detect, ObitTableOTFSolnRow *row)
{
  ofloat gain, diff, sum, w1, w2, delta, d1=0.0, d2=0.0, fblank=ObitMagicF();
  ofloat sum2, sum3, rms, caloff;
  olong i, count;

  /* If no data merely flag solution and return */
  if (nsample<=0) {
    row->add[detect]  = fblank;
    row->mult[detect] = 0.0;
    row->wt[detect]   = 0.0;
    row->cal[detect]  = fblank;
    return;
  } /* end flag if no data */
  
  /* Loop over data */
  count = 0;
  sum = sum2 = sum3 = 0.0;
  for (i=0; i<nsample; i++) {
    if ((cal[i]>0.0) && (data[i]!=fblank)) { /* is this a valid cal value */
      /* First special case */
      if (i==0) {
	if ((cal[i+1]<=0.0) && (data[i+1]!=fblank)) {
	  diff = data[i] - data[i+1];
	  if (diff>0.0) { /* ignore clearly bad values */
	    count++;
	    sum  += calJy / diff;
	    sum2 += (calJy / diff) * (calJy / diff);
	    sum3 += diff;
	  }
	}
	
	/* last datum special case */
      } else if (i==(nsample-1)) {
	if ((cal[i-1]<=0.0) && (data[i-1]!=fblank)) {
	  diff = data[i] - data[i-1];
	  if (diff>0.0) { /* ignore clearly bad values */
	    count++;
	    sum  += calJy / diff;
	    sum2 += (calJy / diff) * (calJy / diff);
	    sum3 += diff;
	  }
	}

	/* neither first nor last */
      } else {
	/* determine interpolation weights */
	delta = time[i+1] - time[i-1]; /* time difference */
	w1 = w2 = 0.0;
	d1 = d2 = 0.0;

	/* prior time */
	if ((cal[i-1]<=0.0) && (data[i-1]!=fblank)) {
	  d1 = data[i] - data[i-1];
 	  if (d1>0.0) { /* ignore clearly bad values */
	    w1 = (time[i+1] - time[i]) / delta;
	  }
	}
	/* following time */
	if ((cal[i+1]<=0.0) && (data[i+1]!=fblank)) {
	  d2 = data[i] - data[i+1];
 	  if (d2>0.0) { /* ignore clearly bad values */
	    w2 = (time[i] - time[i-1]) / delta;
	  }
	}
	
	/* determine gain */
	if (w1+w2 > 0.0) {
	  diff = (w1 * d1 + w2 * d2) / (w1+w2);
	  count++;
	  sum  += calJy / diff;
	  sum2 += (calJy / diff) * (calJy / diff);
	  sum3 += diff;
	}
      } /* end of process by position in time */
    } /* End of cal value */
  } /* end loop over data */
  
  /* determine gain */
  if (count>0) gain = sum / count;
  else gain = fblank;  /* no valid data */
  
  /* determine fractonal RMS */
  if ((count>0) && (gain!=fblank)) {
    sum2 /= count;
    caloff = sum3 / count; /* average cal value */
    rms = sum2 - gain*gain;
    if (rms>0.0) rms = sqrt(rms) / gain;
  } else {
    rms = fblank;  /* no valid data */
    caloff = 0.0;
  }
  
  /* if not acceptable replace with lastgood */
  if (rms > minrms) gain = lastgood;

  /* Set gain in record */
  if ((gain!=fblank) && (gain!=0.0)) gain = 1.0 / gain; /* as correction */
  row->mult[detect] = gain;
  row->cal[detect]  = caloff;
} /* end ObitOTFGetSolnGainSolve  */

/**
 * Return average of 50% of values around median
 * Does sort then returns average of center array
 * \param array array of values, on return will be in ascending order
 * \param incs  increment in data array, >1 => weights
 * \param n     dimension of array
 * \return      Median/average
 */
static ofloat GetSolnMedianAvg (ofloat *array, olong incs, olong n)
{
  olong i, cnt, wid, ngood, ind, hi, lo;
  ofloat temp, wt, out=0.0, fblank = ObitMagicF();
  ofloat maxWt, clipWt, center;
  if (n<=0) return fblank;

  /* If weights, clip to top 95% and set to 2.0e20 */
  ngood = 0;
  if (incs>1) {
    maxWt = -1.0e20;
    for (i=0; i<n; i++) if (array[i*incs]!=fblank) maxWt = MAX(maxWt, array[i*incs+1]);
    /* Clip weights and set to large value */
    clipWt = 0.05*maxWt;
    for (i=0; i<n; i++) {
      if ((array[i*incs]==fblank) || (array[i*incs+1]<clipWt)) {
	array[i*incs]   = 2.0e20;
	array[i*incs+1] = 0.0;
      } else ngood++;
    }
  } else { /* No weights */
    for (i=0; i<n; i++) {
      if (array[i]==fblank) array[i] = 2.0e20;
       else ngood++;
    }
  }
  if (ngood<1) return fblank;  /* Any good data? */

  /* Sort to ascending order */
  qsort ((void*)array, n, MIN(2,incs)*sizeof(ofloat), compare_gfloat);

  /* How many good values? */
  cnt = 0;
  for (i=0; i<n; i++) if (array[i*incs]<1.0e20) cnt++;
  wid = cnt/4;
  ngood = cnt;

  ind = (ngood/incs)/2;  /* make center pixel even to get value not wt */
  ind = (ind/2)*2;
  out = array[ind];

  /* average central 50% values */
  if (ngood>5) {
    center = (ngood-1)/2.0;
    lo = (olong)(center - wid + 0.6);
    hi = (olong)(center + wid);
    /* Weighted? */
    if (incs>1) {
      temp = 0.0; wt = 0.0;
      for (i=lo; i<=hi; i++) {
	if (array[i*incs]<1.0e20) {
	  wt += array[i*incs+1];
	  temp += array[i*incs]*array[i*incs+1];
	}
      }
      if (wt>0.0) out = temp/wt;
    } else { /* unweighted */
      temp = 0.0; cnt = 0;
      for (i=lo; i<=hi; i++) {
	if (array[i]<1.0e20) {
	  cnt++;
	  temp += array[i];
	}
	if (cnt>0) out = temp/cnt;
      }
    } /* end if weighting */
  }

  /* reblank values */
  for (i=0; i<n; i++) {
    if (array[i*incs]>1.0e20) array[i*incs] = fblank;
  }

  return out;
} /* end GetSolnMedianAvg */

/**
 * Fit for instrumental model in each detector
 * offset is the median value, gain determined from cal-on values
 * \param geom   Array geometry structure
 * \param desc   OTF descriptor
 * \param nDet   Number of detectors in data
 * \param nTime  Number of times in data, iCal
 * \param data   OTF Data array, destructive operation
 * \param iCal   If > 0 then the cal value in Jy and this is a cal on.
 * \param row    Solution table row
 */
static void 
ObitOTFGetInstSolve (ObitOTFArrayGeom *geom, ObitOTFDesc *desc, 
		     olong nDet, olong nTime, ofloat **data, ofloat *iCal, 
		     ObitTableOTFSolnRow *row)
{
  olong i, j, count;
  ofloat value, corr, cal=0.0, diff, sum, fblank=ObitMagicF();

  /* Initialize */
  row->dAz = 0.0;  /* no pointing corrections */
  row->dEl = 0.0; /* no pointing corrections */

  /* Get calibration value gain */
  for (i=0; i<nDet; i++) {
    sum = 0.0;
    count = 0;
    for (j=1; j<nTime; j++) {
      if (((iCal[j]>0.0) && (iCal[j-1]<=0.0)) &&
	  ((data[i][j]!=fblank)&&(data[i][j-1]!=fblank))){
	cal = iCal[j];
	diff = data[i][j] - data[i][j-1];
	sum += diff;
	count++;
      }
      if (count>0) {
	row->mult[i] = cal / (sum / count);
	row->cal[i] = sum / count;
      } else {
	row->mult[i] = fblank;
	row->cal[i] = fblank;
      }
    }
  /* Subtract cals */
    if (row->cal[i] != fblank) {
      for (j=0; j<nTime; j++) {
	if (iCal[j]>0.0) data[i][j] -= row->cal[i];
      }
    }
  } /* End loop over detector */

  /* median of data per detector, will reorder data */
  for (i=0; i<nDet; i++) {
    value = GetSolnMedianAvg (data[i], 1, nTime);
    if (value!=fblank) corr = -value;
    else corr = fblank;
    row->add[i]  = corr;
  }


} /* end ObitOTFGetInstSolve */

/**
 * Fit polynomial y = f(poly, x)
 * Use gsl package.
 * \param poly   [out] polynomial coef in order of increasing power of x
 * \param order  order of the polynomial
 * \param x      values at which y is sampled
 * \param y      values to be fitted
 * \param wt     weights for values
 * \param n      number of (x,y) value pairs
 */
static void  FitBLPoly (ofloat *poly, olong order, ofloat *x, ofloat *y, 
			ofloat *wt, olong n)
{
  olong i, j, k, good, p=order+1;
  ofloat fgood=0.0;
  double xi, chisq;
  gsl_matrix *X, *cov;
  gsl_vector *yy, *w, *c;
  gsl_multifit_linear_workspace *work;
  ofloat fblank = ObitMagicF();

  /* Only use good data */
  good = 0;
  for (i=0; i<n; i++) if (y[i]!=fblank) {good++; fgood = y[i];}

  /* If only one good datum - use it */
  if (good==1) {
    poly[0] = fgood;
    poly[1] = 0.0;
  }

  if (good<(order+1)) {  /* Not enough good data */
    poly[0] = fblank;
    return;
  }

  /* If only one and order=0 use it */
  if ((good==1) && (order==0)) {
    poly[0] = y[0];
    return;
  }

  /* allocate arrays */
  X    = gsl_matrix_alloc(good, p);
  yy   = gsl_vector_alloc(good);
  w    = gsl_vector_alloc(good);
  c    = gsl_vector_alloc(p);
  cov  = gsl_matrix_alloc(p, p);
  work = gsl_multifit_linear_alloc (good, p);

  /* set data */
  k = 0;
  for (i=0; i<n; i++) {
    if (y[i]!=fblank) {
      gsl_vector_set(yy, k, y[i]);
      gsl_vector_set(w, k, wt[i]);
      xi = 1.0;
      for (j=0; j<p; j++) {
	gsl_matrix_set(X, k, j, xi);
	xi *= x[i];
      }
      k++;  /* good ones */
    }
  }

  /* Fit */
  gsl_multifit_wlinear (X, w, yy, c, cov, &chisq, work);

  /* get results */
  for (j=0; j<p; j++) poly[j] = gsl_vector_get(c, j);

  /* Deallocate arrays */
  gsl_matrix_free(X);
  gsl_vector_free(yy);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free (work);
} /* end FitBLPoly */

/**
 * Determine median value for each detector (returned as offset) then fit
 * Piecewise linear segments of time length solint.
 * \param solint  desired solution interval in units of x
 * \param npoly   [out] Number of polynomial segments in tpoly, poly
 * \param tpoly   [out] time (x value) of beginning of each segment, 
 *                 tpoly[npoly] = last time in x + epsilon.
 *                 array allocated, must be g_freeed when done
 * \param poly    [out] polynomial (zeroth, first order coef) x npoly
 *                 array allocated, must be g_freeed when done
 * \param offset  [out] offset per detector
 * \param ndetect Number of detectors 
 * \param maxdata Max. number of samples per detector 
 * \param x       values at which y is sampled
 * \param y       values to be fitted
 * \param wt      weights for values
 * \param n       number of (x,y) value pairs
 */
static void  FitMBBLPoly (ofloat solint, olong *npoly, ofloat **tpoly, ofloat **poly, 
			  ofloat *offset, olong ndetect, olong maxdata, 
			  ofloat *x, ofloat *y, ofloat *wt, olong n)
{
  olong i, j, tgood, ngood, np, off, nseg, first;
  ofloat *array=NULL, *tarray=NULL, *wwt=NULL;
  ofloat si, tnext, fblank = ObitMagicF();

  /* Init output */
  *npoly = 0;
  *poly  = NULL;
  *tpoly = NULL;

  /* Median values for each detector */
  tarray = g_malloc0(n*sizeof(ofloat));
  wwt    = g_malloc0(n*sizeof(ofloat));
  for (i=0; i<ndetect; i++) {
    ngood = 0;
    for (j=0; j<n; j++) {
      if ((y[j+i*maxdata]!=fblank) && (wt[j+i*maxdata]>0.0)) {
	tarray[ngood] = y[j+i*maxdata];
	ngood++;
      }
    } /* end accumulate array */
    offset[i] = GetSolnMedianAvg (tarray, 1, ngood);
  } /* End loop over detectors */

  /* Get median residual in each time */
  array = g_malloc0(ndetect*sizeof(ofloat));
  for (j=0; j<n; j++) {
    tgood = 0;
    for (i=0; i<ndetect; i++) {
      if ((y[j+i*maxdata]!=fblank) && (wt[j+i*maxdata]>0.0) && (offset[i]!=fblank)) {
	array[tgood] = y[j+i*maxdata] - offset[i];
	tgood++;
      }
    } /* end accumulate array */
    tarray[j] = GetSolnMedianAvg (array, 1, tgood);
    wwt[j]    = 1.0;
  } /* End loop over detectors */

  /* Create output arrays */
  np = 5.999 + (x[n-1] - x[0]) / solint;  /* How many segments possible? */
  np = MAX (np, 6);
  si = (x[n-1] - x[0]) / (np - 5.0);      /* Make intervals equal */
  *tpoly = g_malloc0(np*sizeof(ofloat));
  *poly  = g_malloc0(2*np*sizeof(ofloat));

  /* Loop over list */
  off   = 0;
  tnext = x[0] + si; /* end of segment */
  nseg  = 0;
  while ((off<(n-1)) && ((nseg+1)<np)) {
    first = off;
    /* Find end of segment */
    while ((x[off+1]<=tnext) && (off<(n-1))) {off++;}
    /* Fit time seqment */
    FitBLPoly (&(*poly)[2*nseg], 1, &x[first], &tarray[first], &wwt[first], off-first);
    /* DEBUG * No common mode
    (*poly)[2*nseg] = 0.0; (*poly)[2*nseg+1] = 0.0; */
    (*tpoly)[nseg] = x[first];
    tnext = x[off+1] + si;  /* end of next segment */
    nseg++;
    /* Ignore orphans */
    if (off>(n-2)) break;
    g_assert (nseg<np); /* Trap overflow */
  } /* end loop over list */

  *npoly = nseg;  /* save number */
  (*tpoly)[nseg] = x[n-1]+0.001*si; /* last plus a bit */

  /* Free internal work arrays */
  g_free(tarray);
  g_free(wwt);
  g_free(array);
} /* end FitMBBLPoly */

/**
 * Remove outlyers from multi-beam data array
 * The Mean and RMS offset from the poly/offset model is computed for 
 * each detectoris computed and then data more discrepant than sigma times
 * the RMS from its detector mean has it's weight set to zero.
 * Halve weight of points with residuals > sigma/2 times the RMS
 * \param npoly   Number of polynomial segments in tpoly, poly
 * \param tpoly   Time (x value) of beginning of each segment, 
 *                tpoly[npoly] = last time+epsilon;
 * \param poly    Polynomial (zeroth, first order coef) x npoly
 * \param offset  offset per detector
 * \param ndetect Number of detectors 
 * \param maxdata Max. number of samples per detector 
 * \param x       values at which y is sampled [n]
 * \param y       values to be fitted [n][ndetect]
 * \param wt      weights for values [n][ndetect]
 * \param n       number of (x,y) value pairs
 * \param sigma   number of sigma from detector mean to flag data
 */
static void  FitMBBLOut (olong npoly, ofloat *tpoly, ofloat *poly, ofloat *offset, 
			 olong ndetect, olong maxdata, 
			 ofloat *x, ofloat *y, ofloat *wt, olong n,
			 ofloat sigma)
{
  olong i, j, off;
  odouble *sum, *count, *mean, *rms, resid, half;
  olong *flagged;
  ofloat atm, fblank = ObitMagicF();

  /* Allocate */
  sum    = g_malloc0(ndetect*sizeof(odouble));
  count  = g_malloc0(ndetect*sizeof(odouble));
  mean   = g_malloc0(ndetect*sizeof(odouble));
  rms    = g_malloc0(ndetect*sizeof(odouble));
  flagged= g_malloc0(ndetect*sizeof(olong));

  /* Get statistics per detector */
  /* Get means */
  for (j=0; j<ndetect; j++) sum[j] = count[j] = 0.0;
  off = 0;
  for (i=0; i<n; i++) {
    /* Find and evaluate time segment */
    while ((tpoly[off+1]<x[i]) && (off<(npoly+1))) {off++;}
    /*atm = poly[2*off] + poly[2*off+1]*(x[i] - tpoly[off]);*/
    atm = poly[2*off] + poly[2*off+1]*(x[i]);
    for (j=0; j<ndetect; j++) {
      if ((y[i+j*maxdata]!=fblank) && (wt[i+j*maxdata]>0.0)) {
	resid = y[i+j*maxdata] - (offset[j] + atm);
	sum[j]   += resid;
	count[j] += 1.0;
      }
    }
  }
  for (j=0; j<ndetect; j++) {
    mean[j] = sum[j]  / MAX (1.0, count[j]);
  }

  /* Get RMSes */
  for (j=0; j<ndetect; j++) sum[j] = count[j] = 0.0;
  off = 0;
  for (i=0; i<n; i++) {
    /* Find and evaluate time segment */
    while ((tpoly[off+1]<x[i]) && (off<(npoly+1))) {off++;}
    /*atm = poly[2*off] + poly[2*off+1]*(x[i] - tpoly[off]);*/
    atm = poly[2*off] + poly[2*off+1]*(x[i]);
    for (j=0; j<ndetect; j++) {
      if ((y[i+j*maxdata]!=fblank) && (wt[i+j*maxdata]>0.0)) {
	resid = y[i+j*maxdata] - (offset[j] + atm) - mean[j];
	sum[j]  += resid*resid;
	count[j] += 1.0;
      }
    }
  }
  for (j=0; j<ndetect; j++) {
    rms[j]  = sqrt (sum[j] / MAX (1.0, count[j]));
    flagged[j] = 0;
  }
 

  /* debug
  j = 0;
  fprintf (stderr,"clip level %d %f %f %f\n", j, sigma*rms[j], rms[j],sigma);
  fprintf (stderr,"   %lf %lf \n", sum[j], count[j]); */

  /* flag outliers */
  half = sigma*0.5;
  off = 0;
  for (i=0; i<n; i++) {
    /* Find and evaluate time segment */
    while ((tpoly[off+1]<x[i]) && (off<(npoly+1))) {off++;}
    /*atm = poly[2*off] + poly[2*off+1]*(x[i] - tpoly[off]);*/
    atm = poly[2*off] + poly[2*off+1]*(x[i]);
    for (j=0; j<ndetect; j++) {
      if ((y[i+j*maxdata]!=fblank) && (wt[i+j*maxdata]>0.0)) {
	resid = y[i+j*maxdata] - (offset[j] + atm);
	if (fabs(resid-mean[j]) > sigma*rms[j]) {
	  wt[i+j*maxdata] = 0.0;
	  flagged[j]++;
	}
	/* Only halve weigh of points more discrepant than sigma/2 */
	if (fabs(resid-mean[j]) > half*rms[j] ) wt[i+j*maxdata] *= 0.5;
	/* debug
	if (j==0) 
	  fprintf (stderr," %d %f %f %lf %f\n",
		   i,y[i+j*maxdata],(offset[j] + atm),resid,wt[i+j*maxdata]); */
      }
    }
  }

  /* DEBUG 
  for (j=0; j<ndetect; j++) {
    if (offset[j]!=fblank) 
	  fprintf (stderr,"DEBUG det %d offset %f mean %f rms %f flagged  %d\n",
		   j, offset[j], mean[j], rms[j], flagged[j]); 
	}   END DEBUG */
  

  /* Deallocate arrays */
  if (sum)   g_free(sum);
  if (count) g_free(count);
  if (mean)  g_free(mean);
  if (rms)   g_free(rms);
  if (flagged) g_free(flagged);
} /* end FitMBBLOut */

/**
 * Plot model and residuals
 * \param npoly   Number of polynomial segments in tpoly, poly
 * \param tpoly   Time (x value) of beginning of each segment, 
 *                tpoly[npoly] = last time+epsilon;
 * \param poly    Polynomial (zeroth, first order coef) x npoly
 * \param offset  offset per detector
 * \param ndetect Number of detectors 
 * \param maxdata Max. number of samples per detector 
 * \param order   order of the polynomial
 * \param x       values at which y is sampled [n]
 * \param y       values to be fitted [n][ndetect]
 * \param wt      weights for values [n][ndetect]
 * \param n       number of (x,y) value pairs
 * \param plotDetect  detector number (1-rel) to plot
 * \param t0      time offset to add to x
 * \param err     Error stack, returns if not empty.
 */
static void PlotMBBL (olong npoly, ofloat *tpoly, ofloat *poly, ofloat *offset, 
		      olong ndetect, olong maxdata, 
		      ofloat *x, ofloat *y, ofloat *wt, olong n,
		      olong plotDetect, ofloat t0, ObitErr *err)
{
  olong i, j, nplot, off=0;
  ofloat *ptime, *pdata, ymax, ymin, yy, atm=0.0, fblank = ObitMagicF();
  ObitPlot *plot = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar strTemp[121];
  gchar *routine = "PlotMBBL";

  /* Asked for something? */
  if ((plotDetect<1) || (plotDetect>=ndetect)) return;

  /* initialize plotting */
  plot = newObitPlot (routine);

  /* Plot labels */
  g_snprintf (strTemp, 120, "MB BL model (line) and data, detector %d", plotDetect);
  dim[0] = strlen (strTemp);
  ObitInfoListAlwaysPut (plot->info, "TITLE", OBIT_string, dim, strTemp);
  strncpy (strTemp, "Time (days)", 120);
  dim[0] = strlen (strTemp);
  ObitInfoListAlwaysPut (plot->info, "XLABEL", OBIT_string, dim, strTemp);
  strncpy (strTemp, "Offset", 120);
  dim[0] = strlen (strTemp);
  ObitInfoListAlwaysPut (plot->info, "YLABEL", OBIT_string, dim, strTemp);

  /* initialize plotting */
  ObitPlotInitPlot (plot, NULL, 15,1,1, err);
  if (err->error) {
    plot = ObitPlotUnref(plot);
    Obit_traceback_msg (err, routine, plot->name);
  }

  /* Allocate arrays */
  ptime   = g_malloc0(n*sizeof(ofloat));
  pdata   = g_malloc0(n*sizeof(ofloat));

  /* Find extrema */
  ymax = -1.0e20;
  ymin =  1.0e20;
  j = plotDetect-1;
  off = 0;
  atm = poly[2*off] + poly[2*off+1]*(x[j]);
  for (i=0; i<n; i++) {
    /*if ((y[i+j*maxdata]!=fblank) && (wt[i+j*maxdata]>0.0)) {*/
    if (y[i+j*maxdata]!=fblank) {
      yy = y[i+j*maxdata] - atm - offset[j];
      ymax = MAX (ymax, y[i+j*maxdata]);
      ymin = MIN (ymin, y[i+j*maxdata]);
    }
  }
  /* Set extrema */
  dim[0] = 1;
  ObitInfoListAlwaysPut(plot->info, "YMAX", OBIT_float, dim, (gpointer*)&ymax);
  ObitInfoListAlwaysPut(plot->info, "YMIN", OBIT_float, dim, (gpointer*)&ymin);

  /* model to plot */
  j = plotDetect-1;
  nplot = 0;
  off = 0;
  for (i=0; i<n; i++) {
    /* Find and evaluate time segment */
    while ((tpoly[off+1]<x[i]) && (off<(npoly+1))) {off++;}
    /*atm = poly[2*off] + poly[2*off+1]*(x[i] - tpoly[off]);*/
    atm = poly[2*off] + poly[2*off+1]*(x[i]);
    ptime[nplot] = x[i]+t0;
    pdata[nplot] = offset[j] + atm;
    nplot++;
  }
  /* plot it */
  ObitPlotXYPlot (plot, -1, nplot, ptime, pdata, err);

  /* get  valid data to plot */
  j = plotDetect-1;
  nplot = 0;
  for (i=0; i<n; i++) {
    if ((y[i+j*maxdata]!=fblank) && (wt[i+j*maxdata]>0.0)) {
      ptime[nplot] = x[i]+t0;
      pdata[nplot] = y[i+j*maxdata];
      nplot++;
    }
  }
  /* plot it */
  ObitPlotXYOver (plot, 2,  nplot, ptime, pdata, err);

  /* get  invalid data to plot */
  j = plotDetect-1;
  nplot = 0;
  for (i=0; i<n; i++) {
    if ((y[i+j*maxdata]!=fblank) && (wt[i+j*maxdata]<=0.0)) {
      ptime[nplot] = x[i]+t0;
      pdata[nplot] = y[i+j*maxdata];
      nplot++;
    }
  }
  /* plot it */
  ObitPlotXYOver (plot, 3,  nplot, ptime, pdata, err);
  ObitPlotFinishPlot (plot, err);
  if (err->error) {
    plot = ObitPlotUnref(plot);
    if (ptime)  g_free(ptime);
    if (pdata)  g_free(pdata);
    Obit_traceback_msg (err, routine, routine);
  }

  /* Deallocate arrays */
  plot = ObitPlotUnref(plot);
  if (ptime)  g_free(ptime);
  if (pdata)  g_free(pdata);
} /* end PlotMBBL */

/**
 * Remove variations in data, calculate cal delta and weight
 * from inverse variance of data.
 * Per detector:
 * 1) Determined cal from averaging each segment 
 *    and differencing transitions;
 * 2) Corrects data for cal value
 * 3) Fit polynomial (5th order) to data
 * 4) Determine rms deviation
 * 5) Clip points > 5 sigma from 0
 * 6) Refit polynomial
 * 7) redetermine rms deviation
 * 8) Flag scans with weights with entries further than 10X low or
          5X high from the median weight
 * \param nDet        Number of detectors
 * \param nTime       Number of times
 * \param iCal        is Cal [nTime] 0 = cal off
 * \param accum       data accumulation[det][time]  modified on return
 * \param time        time per datum[time]  modified on output
 * \param lastCal     [out] Cal signal difference per detector
 * \param lastWeight  [out] Weight per detector (1/variance)
 * \param fitWate     if TRUE then determine weight from this scan, else keep last
 */
static void doPARCalWeight (olong nDet, olong nTime, ofloat *iCal, ofloat **accum, 
			    ofloat *time, ofloat *lastCal, ofloat *lastWeight, 
			    gboolean fitWate)
{
  olong i, j, which, other;
  olong naxis[1], pos[1], msample, order;
  ofloat cal, onCnt, offCnt, fblank = ObitMagicF();
  ofloat cnt[2], sum[2];
  ofloat *tmp = NULL, *wt = NULL, poly[10], rms, t0, medWt;
  ofloat flipCnt, flipSum;
  gboolean isCal, redo, Cal[2]={FALSE,FALSE};
  ObitFArray *tmpFA = NULL;

   /* error checks */
  g_assert(iCal!=NULL);
  g_assert(accum!=NULL);
  g_assert(lastCal!=NULL);
  g_assert(lastWeight!=NULL);

  if ((nDet<=0) || (nTime<=0)) return;  /* Anything to do? */

  /* Work ObitFArray */
  naxis[0] = nTime;
  tmpFA    = ObitFArrayCreate ("temp", 1L, naxis);
  pos[0]   = 0;
  tmp      = ObitFArrayIndex (tmpFA,  pos);  /* Pointer to data array */

  /* Order of polynomial */
  order = MAX(1, nTime/50);
  order = MIN (5, order);

  /* Work Weight array */
  wt = g_malloc0(nTime*sizeof(ofloat));
  for (i=0; i<10; i++) poly[i] = 0.0;

  /* reference times to first */
  t0 = time[0];
  for (i=0; i<nTime; i++) time[i] -= t0;

  /* Loop over detector */
  for (j=0; j<nDet; j++) {
    
    cal = lastCal[j];  /* Last good cal */
    
    /* Difference average On and off */
    onCnt = offCnt =  0.0;
    cnt[0] = cnt[1] = 0.0; sum[0] = sum[1] = 0.0;
    which = 0; other = 1-which;
    flipCnt = flipSum = 0.0;
    for (i=0; i<nTime; i++) {
      if (accum[j][i] != fblank) {
	isCal = iCal[i]==0.0;
	
	/* Did state switch? */
	if (isCal!=Cal[which]) {
	  other = 1-which;  /* not which */
	  
	  /* Have both states */
	  if ((cnt[0]>0.0) && (cnt[1]>0.0)) {
	    flipCnt++;
	    /* Which way did it go? */
	    if (isCal) { /* off to on */
	      flipSum += (sum[which]/cnt[which]) - (sum[other]/cnt[other]);
	    } else {     /* on to off */
	      flipSum += (sum[other]/cnt[other]) - (sum[which]/cnt[which]);
	    }
	  }
	  
	  /* Flip which */
	  which = other;
	  cnt[which] = 0.0;
	  sum[which] = 0.0;
	  Cal[which] = isCal;
	}

	/* Sums for this cycle part */
	cnt[which]++;
	sum[which] += accum[j][i];
	/* Total counts */
	if (isCal) {  /* cal on */
	  onCnt++;
	} else {   /* cal off */
	  offCnt++;
	}
      }
    } /* end loop over time averaging cal */
    
    /* If no valid data skip this detctor */
    if ((onCnt<1) && (offCnt<1)) continue;
    
    if (flipCnt>0) { /* Need both cal on and cal off */
      cal = flipSum/flipCnt;
      lastCal[j] = cal; 

      /* Correct cal on to cal off */
      for (i=0; i<nTime; i++) {
	if ((accum[j][i] != fblank) && (iCal[i]!=0.0)) accum[j][i] -= cal;
      } /* end loop over time correcting cal */
    }  /* End both cal on and cal off */

    /* Set fitting weights */
    for (i=0; i<nTime; i++) {
      if (accum[j][i] != fblank) wt[i] = 1.0;
      else wt[i] = 0.0;
    }

    /* Need to determine weight? */
    if (!fitWate) continue;
      
    /* Flatten curve - first fit 5th order polynomial */
    msample = nTime;
    FitBLPoly (poly, order, time, accum[j], wt, msample); 
    
    /* determine RMS from histogram analysis - use FArray class */
    for (i=0; i<nTime; i++) {
      if (accum[j][i]!=fblank) 
	tmp[i] = accum[j][i] - EvalPoly (order, poly, time[i]);
      else tmp[i] = fblank;
    }
    rms =ObitFArrayRawRMS (tmpFA);
    if ((rms!=fblank) && (rms>0.0)) {
      rms = MAX (rms, 1.0);
      lastWeight[j] = 1.0e-4*cal*cal / (rms*rms);
    }

    /* DEBUG
       if ((lastWeight[j]!=fblank) && (lastWeight[j]>20.0)) {
       fprintf (stderr,"DEBUG detector %d weight %f rms %f cal %f\n",
       j, lastWeight[j],rms,cal);
    } */

    /* Clip >5 sigma deviations */
    redo = FALSE;
    for (i=0; i<nTime; i++) {
      if ((tmp[i]!=fblank) && (fabs(tmp[i]) > 5.0*rms)) {
	wt[i] = 0.0;
	tmp[i] = fblank;
	redo = TRUE;
      }
    }

    if (redo) {
      /* Refit curve */
      FitBLPoly (poly, order, time, accum[j], wt, msample); 
      
      /* Redetermine RMS */
      for (i=0; i<nTime; i++) {
	if (accum[j][i]!=fblank) 
	  tmp[i] = accum[j][i] - EvalPoly (order, poly, time[i]);
	else tmp[i] = fblank;
      }
      rms = ObitFArrayRawRMS (tmpFA);
      if ((rms!=fblank) && (rms>0.0)) {
	rms = MAX (rms, 1.0);
	lastWeight[j] =  1.0e-4*cal*cal / (rms*rms);
      }

    }/* End redo after clipping */
  } /* end loop over detectors */

  /* flag detectors with cal amplitude more than 10 X from median */
  if (wt) g_free(wt);
  wt = g_malloc0(nDet*2*sizeof(ofloat));
  for (j=0; j<nDet; j++) {
    if (lastCal[j]!=fblank) wt[j] = fabs(lastCal[j]);
    else wt[j] = fblank;
  }
  medWt = medianValue (wt, 1, (olong)nDet);
  if (medWt!=fblank) {
    for (j=0; j<nDet; j++) {
      if ((fabs(lastCal[j])<0.1*medWt) || (fabs(lastCal[j])>10.0*medWt)) {
	/* Bad */
	/* lastCal will blank lastWeight[j] = 0.0;*/
	lastCal[j]    = fblank;
      }
    }
  }

  /* flag  weights more than 10 X low or 5 X high from median */
  if (fitWate){
    for (j=0; j<nDet; j++) wt[j] = lastWeight[j];
    medWt = medianValue (wt, 1, (olong)nDet);
    if (medWt!=fblank) {
      for (j=0; j<nDet; j++) {
	if ((lastWeight[j]<0.1*medWt) || (lastWeight[j]>5.0*medWt)) {
	  /* Bad */
	  lastWeight[j] = 0.0;
	}
      }
    }
  }

  /* Cleanup */
  if (wt) g_free(wt);
  ObitFArrayUnref(tmpFA);
} /* end doPARCalWeight */

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

/**
 * Determine the zenith opacity from either a constant scalar or
 * an array of (time,tau0) values.
 * Values MUST be in increasing time order.
 * \param taudim   Dimensionality of taudata
 * \param taudata  scalar tau0 or array of (time,tau0)
 * \return Zenith opacity
 */
static ofloat getTau0 (ofloat time, gint32 taudim[MAXINFOELEMDIM], ofloat *taudata)
{
  ofloat out = 0.0;
  ofloat wt1, wt2;
  olong i, ntime, pre=0;

  if (taudata==NULL) return out; /* Given anything? */

  /* Scalar?*/
  if (taudim[0]==1)  out = *taudata;
  else if (taudim[0]==2) {
    if (taudim[1]==1) out = taudata[1];  /* Only one entry */
    else {  /* Multiple entries */
      ntime = taudim[1];
      /* Interpolate/extrapolate */
      if (time<=taudata[0])  out = taudata[1];  /* Before first */
      if (time>=taudata[(ntime-1)*2])  out = taudata[2*(ntime-1)+1];  /* After last */
      else {
	/* interpolate */
	/* find last preceeding */
	for (i=0; i<ntime; i++) if (time>taudata[2*i]) pre = i;
	wt2 = (time - taudata[2*pre]) / (taudata[2*ntime-2] - taudata[0]);
	wt1 = 1.0 - wt2;
	out = wt1*taudata[2*pre+1] + wt2*taudata[2*pre+3];
      }
    }
  }
  return out;
} /* end getTau0*/

/**
 * Get RMS values of an array.
 * \param array array of values may be magic value blanked
 * \param n     dimension of array
 * \return      RMS about mean value, possibly blanked
 */
static ofloat RMSValue (ofloat *array, olong n)
{
  olong i, count;
  ofloat sum, sum2, rawRMS, fblank = ObitMagicF();

  count = 0; sum = sum2 = 0.0; 
  for (i=0; i<n; i++) if (array[i]!=fblank) {
    count++;
    sum  += array[i];
    sum2 += array[i] * array[i];
  }
  if (count<5) return fblank; /* better have some */

  /* Get raw RMS */
  rawRMS = (sum2/count) - ((sum / count) * (sum / count));
  if (rawRMS>0.0) rawRMS = sqrt(rawRMS);
  return rawRMS;
} /* end RMSValue */

/**
 * Lookup the ID number of Target named "Blank"
 * \param inOTF    Input OTF data
 * \param err      Error stack, returns if not empty.
 * \return "Blank" Target ID if found, else -1.0
 */
static ofloat FindBlankID (ObitOTF *inOTF, ObitErr *err)
{
  ofloat out = -1.0;
  ObitTableOTFTarget *TarTable=NULL;
  gint32 dim[2];
  olong iver, target[2];
  gchar *Blank = "Blank";
  gchar *routine = "FindBlankID";
 

  if (err->error) return out; /* existing error? */

     iver = 1;
    TarTable = newObitTableOTFTargetValue (inOTF->name, (ObitData*)inOTF, &iver, 
					   OBIT_IO_ReadOnly, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, out);
    if (TarTable==NULL) {/*No Target table */
      out = -1.0;
    } else {  /* Table exists */
      dim[0] = strlen(Blank);
      dim[1] = 1;
       /* Do lookup */
      ObitTableOTFTargetLookup (TarTable, dim, Blank, target, err);
      if(err->error)  Obit_traceback_val (err, routine, inOTF->name, out);
      TarTable = ObitTableOTFTargetUnref(TarTable); /* release table */
      out = (ofloat)target[0];
   }
 return out;
} /* end FindBlankID */
