/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2013                                          */
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
#include "ObitOTFCalUtil.h"
#include "ObitOTFGetSoln.h"
#include "ObitTableOTFIndex.h"
#include "ObitTableOTFTargetUtil.h"
#include "ObitTableOTFFlag.h"
#include "ObitTableOTFBP.h"
#include "ObitPlot.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_multifit_nlin.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFCalUtil.c
 * ObitOTF class calibration utility function definitions.
 */

/*---------------Private function prototypes----------------*/
/** Private: Calibrator fit for single scan */
static void FitCalScan (ObitOTF *inOTF, olong detect, olong scan,
			ofloat *TRx, ofloat *CalJy, 
			ofloat *RAoff, ofloat *Decoff, ofloat *TimeOff,
			ofloat *data, ofloat *elev, ObitErr *err);

/** Private: Average values for a single scan */
static void FitCalAverage (ObitOTF *inOTF, olong detect, olong scan,
			ofloat *calOff, ofloat *cal, 
			ofloat *data, olong *targID, ObitErr *err);

/** Private: Average spectra for a single scan */
static void SpectrumAverage (ObitOTF *inOTF, olong scan,  ofloat *data, 
			     ofloat *work, olong *targID, ofloat *parAng, 
			     ofloat *Time, ObitErr *err);

//** Private: Average onsource and offsource for a single scan */
static void FitCalNod (ObitOTF *inOTF, olong detect, olong scan,
		       ofloat *avgOff, ofloat *avgOn, ofloat *avgCal, 
		       olong *targID, ObitErr *err);

/** Private: Convert list of positions and times to those wrt calibrator */
static void SetWrtCal (odouble RACal, odouble DecCal, 
		       olong Count, ofloat *RA, ofloat *Dec, ofloat *Time);

/** Private: Fit Gaussian to calibrator scan */
static void FitCalGauss (olong Count, ofloat *Time, ofloat *Temp, 
			 ofloat *base, ofloat *peak, ofloat *center, 
			 ofloat *sigma, ObitErr *err);

/** Private: Read and fit tipping scan */
static void FitTipScan (ObitOTF *inOTF, olong scan, ofloat minEl, 
			ofloat *TRx, ofloat *ATemp, ofloat *RMS, ObitErr *err);

/** Private: fit tipping data */
static ofloat FitTip (olong ndata, ofloat *el, ofloat *Temp, ofloat *Tatm, 
		      ofloat *Trx);

/** Private: plot cal scan */
static void PlotCalGauss (olong Count, ofloat *Time,
			  olong scan, olong detect, ofloat *Temp,
			  ofloat base, ofloat peak, ofloat center, 
			  ofloat sigma, ObitErr *err);

/** Private: plot tip scan */
static void PlotTip (olong Count, olong scan, olong detect, ofloat tau0,
		     ofloat *airmass, ofloat *Temp, ofloat Tatm,
		     ofloat Trx, ObitErr *err);


#define MAXSAMPLE 5000   /* Maximum number of samples in a calibrator scan */
#define MAXTIP   10000   /* Maximum number of samples in a tipping scan */
#define ATEMP    290.0   /* Default effective atmospheric Temp */

/*----------------------Public functions---------------------------*/
/**
 * Fits calibrator scan.
 * This routines presumes that the selected detector(s) is scanned 
 * linearly through the expected position of the calibrator .
 * Gets calibrator information from the Target table (position, flux density)
 * Input values on inOTF
 * \li "Scan"     OBIT_Int (?,1,1) Scan numbers to process
 * \li "Tau0"     OBIT_float (1,1,1) Zenith opacity in nepers [def 0].
 * \li "aTemp"    OBIT_float (*,1,1) Effective atmospheric temperature in cal units.
 *                i.e. the additional offset per airmass due to the atmosphere.
 *                per detector in units of the cal. [def 0]
 * \li "doPlot"   OBIT_Bool (1,1,1) If present and true, plot data and models [def FALSE]
 * Output values on inOTF
 * \li "tRx"      OBIT_float (*,1,1) Receiver temperature per detector in units of the cal
 *                -1 => not determined
 * \li "calJy"    OBIT_float (*,1,1) Noise cal value in Jy, per detector
 *                -1 => not determined
 * \li "RAoff"    OBIT_float (*,1,1) Offset in deg to add to RA, 
 *                Not implemented as the correction applied by M&C
 * \li "Decoff"   OBIT_float (*,1,1) Offset in deg to add to Dec
 *                Not implemented as the correction applied by M&C
 * \li "Timeoff"  OBIT_float (*,1,1) Offset in time (day) (actual-expected peak)
 * \param inOTF    Input OTF data. Maximum MAXSAMPLE integrations per scan.
 * \param detect   detector number (0-rel), -1 => all
 * \param err      Error stack, returns if not empty.
 */
void ObitOTFCalUtilFitCal (ObitOTF *inOTF, olong detect, ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat tau0, ftemp, *farr=NULL, *ATemp=NULL, *calJy=NULL, *TRx=NULL;
  ofloat *RAoff=NULL, *Decoff=NULL, *Timeoff=NULL, fblank = ObitMagicF();
  olong i, cnt, ndetect,  mdetect, loDet, hiDet;
  olong scans[4], iscan, nscan, scan;
  ofloat *data=NULL, sTRx, scalJy, sTimeoff, airmass, elev;

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));

  /* How many detectors? */
  ndetect = inOTF->geom->numberDetect;
  /* How many to use */
  if (detect<0) {
    mdetect = ndetect;
    loDet = 0;
    hiDet = ndetect-1;
  } else {
    mdetect = 1;
    loDet = detect;
    hiDet = detect;
  }

  /* Create arrays */
  farr    = g_malloc0(ndetect*sizeof(float));
  data    = g_malloc0(MAXSAMPLE*ndetect*sizeof(float));
  ATemp   = g_malloc0(ndetect*sizeof(float));
  calJy   = g_malloc0(4*ndetect*sizeof(float));
  TRx     = g_malloc0(4*ndetect*sizeof(float));
  RAoff   = g_malloc0(4*ndetect*sizeof(float));
  Decoff  = g_malloc0(4*ndetect*sizeof(float));
  Timeoff = g_malloc0(4*ndetect*sizeof(float));

  /* Get control parameters */
  /* scan - no default */
  ObitInfoListGet (inOTF->info, "Scan", &type, dim,  (gpointer)&scans, err);
  nscan = dim[0];
 
  /* Default tau0 = 0.0 */
  ftemp = 0.0; type = OBIT_float; dim[0] = 1;
  ObitInfoListGetTest(inOTF->info, "Tau0", &type, dim, (gpointer*)&ftemp);
  tau0 = ftemp;

  /* Default ATemp = 0 */  
  for (i=0; i<ndetect; i++) farr[i] = 0.0;
  type = OBIT_float; dim[0] = ndetect;
  ObitInfoListGetTest(inOTF->info, "Tau0", &type, dim, (gpointer*)&ftemp);
  for (i=0; i<ndetect; i++) ATemp[i] = farr[i];

  /* Loop over input scans */
  for (iscan = 0; iscan<nscan; iscan++) {
    scan = scans[iscan];
    FitCalScan (inOTF, detect, scan, &TRx[iscan*ndetect], &calJy[iscan*ndetect], 
		&RAoff[iscan*ndetect], &Decoff[iscan*ndetect], &Timeoff[iscan*ndetect], 
		data, &elev, err);

    /* Make corrections for atmosphere */
    /* Number of atmospheres seen through */
    if (elev > 0.1) 
      airmass = 1.0 / cos (DG2RAD * (90.0 - elev));
    else
      airmass = 10.0;
    for (i=loDet; i<=hiDet; i++) {
      calJy[iscan*ndetect+i] *= exp (-tau0 * airmass);
      TRx[iscan*ndetect+i]   -= ATemp[i] * airmass;
    }
    
  } /* end loop over scans */
  
  /* Average TRx, calJy */
  for (i=loDet; i<=hiDet; i++) {
    sTRx = scalJy = sTimeoff = 0.0;
    cnt = 0;
    for (iscan = 0; iscan<nscan; iscan++) {
      if (TRx[iscan*ndetect+i]!=fblank) {
	sTRx     += TRx[iscan*ndetect+i];
	scalJy   += calJy[iscan*ndetect+i];
	sTimeoff += Timeoff[iscan*ndetect+i];
	cnt++;
      }
    }
    if (cnt>0) {
      TRx[i]     = sTRx / cnt;
      calJy[i]   = scalJy / cnt;
      Timeoff[i] = sTimeoff / cnt;
    } else {
      TRx[i]     = -1.0;
      calJy[i]   = -1,0;
      Timeoff[i] =  0.0;
    }
  } /* end average over scans */

  /* Save values on inOTF */
  dim[0] = ndetect; dim[1] = 1;
  ObitInfoListPut (inOTF->info, "tRx",     OBIT_float, dim,  (gpointer)TRx, err);
  ObitInfoListPut (inOTF->info, "calJy",   OBIT_float, dim,  (gpointer)calJy, err);
  ObitInfoListPut (inOTF->info, "RAoff",   OBIT_float, dim,  (gpointer)RAoff, err);
  ObitInfoListPut (inOTF->info, "Decoff",  OBIT_float, dim,  (gpointer)Decoff, err);
  ObitInfoListPut (inOTF->info, "Timeoff", OBIT_float, dim,  (gpointer)Timeoff, err);

  /* Cleanup */
  if (farr)    g_free(farr);
  if (data)    g_free(data);
  if (ATemp)   g_free(ATemp);
  if (calJy)   g_free(calJy);
  if (TRx)     g_free(TRx);
  if (RAoff)   g_free(RAoff);
  if (Decoff)  g_free(Decoff);
  if (Timeoff) g_free(Timeoff);
 
}  /* end ObitOTFCalUtilFitCal */

/**
 * Fits calibrator On/Off scan pair.
 * Gets calibrator information from the Target table (flux density)
 * Input values on inOTF
 * \li "Scan"     OBIT_Int (2,1,1) Scan numbers to process
 * Output values on inOTF
 * \li "tRx"      OBIT_float (*,1,1) Receiver temperature per detector in units of the cal
 *                average "off" in units of the cal.
 * \li "calJy"    OBIT_float (*,1,1) Noise cal value in Jy, per detector
 * \param inOTF    Input OTF data. Maximum MAXSAMPLE integrations per scan.
 * \param detect   detector number (0-rel), -1 => all
 * \param err      Error stack, returns if not empty.
 */
void ObitOTFCalUtilFitOnOff (ObitOTF *inOTF, olong detect, ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ObitTableOTFTarget* targetTable=NULL;
  ofloat *farr=NULL, *calJy=NULL, *TRx=NULL;
  ofloat *calOff=NULL, *calAvg=NULL;
  olong i, ndetect,  mdetect, loDet, hiDet;
  olong scans[2], iscan, nscan, scan, targID[2]={0,0};
  olong ver;
  odouble RACal, DecCal;
  ofloat FluxCal, fblank = ObitMagicF();
  ofloat *data=NULL, delta, aveCal;
  gchar *routine = "ObitOTFCalUtilFitOnOff";

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));

  /* How many detectors? */
  ndetect = inOTF->geom->numberDetect;
  /* How many to use */
  if (detect<0) {
    mdetect = ndetect;
    loDet = 0;
    hiDet = ndetect-1;
  } else {
    mdetect = 1;
    loDet = detect;
    hiDet = detect;
  }

  /* Create arrays */
  farr    = g_malloc0(ndetect*sizeof(float));
  data    = g_malloc0(MAXSAMPLE*ndetect*sizeof(float));
  calOff  = g_malloc0(2*ndetect*sizeof(float));
  calAvg  = g_malloc0(2*ndetect*sizeof(float));
  calJy   = g_malloc0(ndetect*sizeof(float));
  TRx     = g_malloc0(ndetect*sizeof(float));

  /* Get control parameters */
  /* scan - no default */
  ObitInfoListGet (inOTF->info, "Scan", &type, dim,  (gpointer)&scans, err);
  nscan = MIN (2, dim[0]);
 
  /* Loop over input scans */
  for (iscan = 0; iscan<nscan; iscan++) {
    scan = scans[iscan];
    FitCalAverage (inOTF, detect, scan, &calOff[iscan*ndetect], &calAvg[iscan*ndetect], 
		data, &targID[iscan], err);
    if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  } /* end loop over scans */

  /* Get cal flux density */
  ver = 1;
  targetTable = 
    newObitTableOTFTargetValue ("TargetTable", (ObitData*)inOTF, &ver, OBIT_IO_ReadWrite, 
				err);
  ObitTableOTFTargetGetSource (targetTable, targID[0], &RACal, &DecCal, &FluxCal, err);
  targetTable = ObitTableOTFTargetUnref(targetTable);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Make sure there is something */
  if (FluxCal==0.0) {
    Obit_log_error(err, OBIT_Error, "%s: MISSING Calibrator info for %d %f in %s", 
		   routine, targID[0], FluxCal, inOTF->name);
  }
  
  /* Average TRx, calJy */
  for (i=loDet; i<=hiDet; i++) {
    if ((calOff[i]!=fblank) && (calOff[i+ndetect]!=fblank)) {
      delta  = calOff[i] - calOff[i+ndetect];
    } else {
      delta  = fblank;
    }
    if ((calAvg[i]!=fblank) && (calAvg[i+ndetect]!=fblank)) {
      aveCal = 0.5 * (calAvg[i] + calAvg[i+ndetect]);
    } else {
      aveCal  = fblank;
    }
    if ((delta!=fblank) && (aveCal!=fblank)) {
      TRx[i]     = calOff[i+ndetect] / aveCal;
      calJy[i]   = FluxCal / (delta / aveCal);
    } else {
      TRx[i]     = fblank;
      calJy[i]   = fblank;
    }
  } /* end loop over detectors */

  /* Save values on inOTF */
  dim[0] = ndetect; dim[1] = 1;
  ObitInfoListPut (inOTF->info, "tRx",     OBIT_float, dim,  (gpointer)TRx, err);
  ObitInfoListPut (inOTF->info, "calJy",   OBIT_float, dim,  (gpointer)calJy, err);

  /* Cleanup */
  if (farr)    g_free(farr);
  if (data)    g_free(data);
  if (calJy)   g_free(calJy);
  if (TRx)     g_free(TRx);
  if (calOff) g_free(calOff);
  if (calAvg) g_free(calAvg);
 
}  /* end ObitOTFCalUtilFitOnOff */

/**
 * Fits calibrator bandpass from On/Off scan pair.
 * Gets calibrator information from the Target table (flux density)
 * Writes BP table
 * \param inOTF    Input OTF data. 
 *    Optional parameters on info:
 * \li "calFlux"    OBIT_float scalar Flux density (Jy) at reference freq
 *                  Default is to get from the Target table
 * \li "calIndex"  OBIT_float scalar Calibrator spectral index def [0]
 * \li "calFPol"   OBIT_float scalar Calibrator poln. orientation in deg [def 0]
 * \li "calEVPA"   OBIT_float scalar Calibrator fractional polarization [def 0]
 * \li "calRM"     OBIT_float scalar Calibrator RM in Rad/m*2 [def 0]
 * \param offScan  off source scan
 * \param onScan   on source scan
 * \param BPVer    BP table version number, 0=>new
 * \param err      Error stack, returns if not empty.
 */
void ObitOTFCalUtilFitBPOnOff (ObitOTF *inOTF,  olong offScan, olong onScan, 
			       olong BPVer, ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ObitTableOTFTarget* targetTable=NULL;
  ObitTableOTFBP* BPTable=NULL;
  ObitTableOTFBPRow* BPRow=NULL;
  ObitOTFDesc *desc = inOTF->myDesc;
  ofloat *data[2]={NULL, NULL}, *wt[2]={NULL,NULL}, Time[2];
  ofloat calFlux=0.0, calIndex=0.0, calFPol=0.0, calEVPA=0.0, calRM=0.0;
  ofloat QPol, UPol, polAdd=0.0, EVPA, lambda, lambda0, deltaNu, refPixNu;
  olong ndetect, iDet, nfeed, nstok, nchan, ifeed, istok, ichan, iRow;
  olong scans[2], sscans[2], iscan, scan, targID[2]={0,0};
  olong ver;
  gboolean isLinear;
  odouble RACal, DecCal, refFreq, freq;
  ofloat FluxCal, freqFlux, parAng=0.0, fblank = ObitMagicF();
  gchar *routine = "ObitOTFCalUtilFitOnOff";

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));

  /* How many detectors? */
  nfeed = desc->inaxes[desc->jlocfeed];
  nstok = desc->inaxes[desc->jlocs];
  /* Stokes limited to 2 solutions */
  if (nstok>2) nstok = 2;
  nchan = desc->inaxes[desc->jlocf];
  ndetect = nfeed * nstok * nchan;

  /* Create arrays */
  data[0] = g_malloc0(ndetect*sizeof(float));
  data[1] = g_malloc0(ndetect*sizeof(float));
  wt[0]   = g_malloc0(ndetect*sizeof(float));
  wt[1]   = g_malloc0(ndetect*sizeof(float));
 
  /* Loop over input scans averaging spectra*/
  sscans[0] = offScan; sscans[1] = onScan;
  for (iscan = 0; iscan<2; iscan++) {
    scan = sscans[iscan];
    /* Set selected scan */
    scans[0] = scan; scans[1] = scan;
    dim[0] = 2;
    ObitInfoListAlwaysPut(inOTF->info, "Scans", OBIT_long, dim, scans);

    SpectrumAverage (inOTF, scan, data[iscan], wt[iscan], &targID[iscan], 
		     &parAng, &Time[iscan], err);
    if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
  } /* end loop over scans */

  /* Check that same source */
  Obit_return_if_fail ((targID[1]==targID[0]), err,
		       "s:Target source not the same for on (%d), off(%d)",
		       targID[1], targID[0]);

  /* Get cal flux density/spectral index */
  ObitInfoListGetTest(inOTF->info, "calIndex", &type, dim, &calIndex);
  ObitInfoListGetTest(inOTF->info, "calFlux",  &type, dim, &calFlux);
  ObitInfoListGetTest(inOTF->info, "calFPol",  &type, dim, &calFPol);
  ObitInfoListGetTest(inOTF->info, "calEVPA",  &type, dim, &calEVPA);
  ObitInfoListGetTest(inOTF->info, "calRM",    &type, dim, &calRM);
  if (calFlux>0.0) {
    FluxCal = calFlux;
  } else { 
    /* Not provided, Get cal flux density */
    ver = 1;
    targetTable = 
      newObitTableOTFTargetValue ("TargetTable", (ObitData*)inOTF, &ver, 
				  OBIT_IO_ReadWrite, err);
    ObitTableOTFTargetGetSource (targetTable, targID[0], &RACal, &DecCal, 
				 &FluxCal, err);
    targetTable = ObitTableOTFTargetUnref(targetTable);
    if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
  } /* End lookup from table */
    /* Make sure there is something */
  if (FluxCal==0.0) {
    Obit_log_error(err, OBIT_Error, "%s: MISSING Calibrator info for %d %f in %s", 
		   routine, targID[0], FluxCal, inOTF->name);
  }

  /* Frequency info */
  refFreq  = desc->crval[desc->jlocf];
  deltaNu  = desc->cdelt[desc->jlocf];
  refPixNu = desc->crpix[desc->jlocf];
  /* Circular or linear feeds */
  isLinear = desc->crval[desc->jlocs]<-4.0;

  /* Calculate the calibration (to Jy) for each channel, leave in data[0] */
  iDet = 0;
  /* Loop over feed */
  for (ifeed=0; ifeed<nfeed; ifeed++) {
    /* Loop over Stokes */
    for (istok=0; istok<nstok; istok++) {
      /* Loop over Channel */
      for (ichan=0; ichan<nchan; ichan++) {
	/* Get frequency/poln related */
	freq     = refFreq + (ichan + 1.0 - refPixNu) * deltaNu;
	freqFlux = FluxCal * pow((freq/refFreq), calIndex);

	/* Corrections for linear feeds */
	if (isLinear) {
	  lambda0 = VELIGHT/refFreq;   /* reference wavelength */
	  lambda  = VELIGHT/freq;      /* channel wavelength */
	  /* Correct EVPA for RM */
	  EVPA = calEVPA*DG2RAD + (lambda*lambda-lambda0*lambda0)*calRM;
	  /* Polarization additions including parallactic angle */
	  QPol = freqFlux * calFPol * cos(2.0*EVPA) * cos(2*parAng*DG2RAD);
	  UPol = freqFlux * calFPol * sin(2.0*EVPA) * sin(2*parAng*DG2RAD);
	  if (istok==0) { /* XX */
	    polAdd = QPol + UPol;
	  } else {        /* YY */
	    polAdd = -QPol - UPol;
	  }
	} /* End isLinear */
	if ((data[0][iDet]!=fblank) && (data[1][iDet]!=fblank) &&
	    (fabs(data[1][iDet]-data[0][iDet])>fabs(0.01*(data[1][iDet])))) {
	  data[0][iDet] = (freqFlux+polAdd) / (data[1][iDet]-data[0][iDet]);
	} else {  /* Bad */
	  data[0][iDet] = fblank;
	}
	iDet++;
      } /* end channel loop */
    } /* end Stokes loop */
  } /* end feed loop */

  /* Write in OTFBP table */
  ver = BPVer;
  if (ver<=0) ver = ObitTableListGetHigh (inOTF->tableList, "OTFBP");
  BPTable = newObitTableOTFBPValue("BP Table", (ObitData*)inOTF, &ver, 
				   OBIT_IO_ReadWrite, nchan, nstok, nfeed, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Open table */
  ObitTableOTFBPOpen (BPTable, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Create Table Row */
  BPRow = newObitTableOTFBPRow (BPTable);
  
  /* Attach  row to output buffer */
  ObitTableOTFBPSetRow (BPTable, BPRow, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* If there are entries in the table, mark it unsorted */
  if (BPTable->myDesc->nrow>0) 
    {BPTable->myDesc->sort[0]=0; BPTable->myDesc->sort[1]=0;}
  
  /* Fill in BP row */
  BPRow->Time   = 0.5*(Time[1]+Time[0]);
  BPRow->TimeI  = fabs(Time[1]-Time[0]);
  BPRow->Target = targID[1];

  for (iDet=0; iDet<ndetect; iDet++) {
    BPRow->mult[iDet] = data[0][iDet];
    BPRow->wt[iDet]   = wt[0][iDet]+wt[1][iDet];
    }
  
  /* write row */
  iRow = BPTable->myDesc->nrow+1;
  ObitTableOTFBPWriteRow (BPTable, iRow, BPRow, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
   
  /* Close BP table */
  ObitTableOTFBPClose (BPTable, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Cleanup */
  BPTable = ObitTableOTFBPUnref(BPTable);
  BPRow   = ObitTableOTFBPRowUnref(BPRow);
  if (data[0])    g_free(data[0]);
  if (data[1])    g_free(data[1]);
  if (wt[0])      g_free(wt[0]);
  if (wt[1])      g_free(wt[1]);
 
}  /* end ObitOTFCalUtilFitBPOnOff */

/**
 * Fits calibrator Nodding scan.
 * Gets calibrator flux density information either from info item
 * \li "calFlux"  OBIT_float (1,1,1) Calibrator flux density at ref. freq
 * \li "calIndex" OBIT_float (1,1,1) Spectral index, [default -0.7]
 * or from the Target table
 * Input values on inOTF
 * \li "Scan"     OBIT_int (1,1,1) Scan number to process
 * Output values on inOTF
 * \li "tRx"      OBIT_float (*,1,1) Receiver temperature per detector in units of the cal
 *                average "off" in units of the cal.
 * \li "calJy"    OBIT_float (*,1,1) Noise cal value in Jy, per detector
 * \param inOTF    Input OTF data. Maximum MAXSAMPLE integrations per scan.
 * \param detect   detector number (0-rel), -1 => all
 * \param err      Error stack, returns if not empty.
 */
void ObitOTFCalUtilFitNod (ObitOTF *inOTF, olong detect, ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ObitTableOTFTarget* targetTable=NULL;
  ofloat *farr=NULL, *calJy=NULL, *TRx=NULL;
  ofloat *avgOff=NULL, *avgOn=NULL, *avgCal=NULL;
  ofloat calFlux, calIndex, chFlux;
  olong i, ndetect,  mdetect, loDet, hiDet, ichan, nchan;
  olong scan, targID=0;
  olong ver;
  odouble RACal, DecCal, refFreq, freqRatio;
  ofloat FluxCal, fblank = ObitMagicF();
  ofloat delta;
  gchar *routine = "ObitOTFCalUtilFitNod";

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));

  /* How many detectors? */
  ndetect = inOTF->geom->numberDetect;
  /* How many to use */
  if (detect<0) {
    mdetect = ndetect;
    loDet = 0;
    hiDet = ndetect-1;
  } else {
    mdetect = 1;
    loDet = detect;
    hiDet = detect;
  }

  /* Create arrays */
  farr    = g_malloc0(ndetect*sizeof(float));
  avgOff  = g_malloc0(ndetect*sizeof(float));
  avgOn   = g_malloc0(ndetect*sizeof(float));
  avgCal  = g_malloc0(ndetect*sizeof(float));
  calJy   = g_malloc0(ndetect*sizeof(float));
  TRx     = g_malloc0(ndetect*sizeof(float));

  /* Get control parameters */
  /* scan - no default */
  ObitInfoListGet (inOTF->info, "Scan", &type, dim,  &scan, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
 
  /* Average scan */
  FitCalNod (inOTF, detect, scan, avgOff, avgOn, avgCal, &targID, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
  
  /* Get cal flux density/spectral index */
  calIndex = -0.7; 
  ObitInfoListGetTest(inOTF->info, "calIndex", &type, dim, &calIndex);
  calFlux = -1.0; 
  ObitInfoListGetTest(inOTF->info, "calFlux", &type, dim, &calFlux);
  if (calFlux>0.0) {
    FluxCal = calFlux;
  } else { 
    /* lookup in target table */
    ver = 1;
    targetTable = 
      newObitTableOTFTargetValue ("TargetTable", (ObitData*)inOTF, &ver, OBIT_IO_ReadWrite, 
				  err);
    ObitTableOTFTargetGetSource (targetTable, targID, &RACal, &DecCal, &FluxCal, err);
    targetTable = ObitTableOTFTargetUnref(targetTable);
    if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
  } /* end lookup in table */

  /* Make sure there is something */
  if (FluxCal==0.0) {
    Obit_log_error(err, OBIT_Error, "%s: MISSING Calibrator info for %d %f in %s", 
		   routine, targID, FluxCal, inOTF->name);
  }
  
  /* reference frequency */
  refFreq = (inOTF->myDesc->crval[inOTF->myDesc->jlocf] + 
	     inOTF->myDesc->cdelt[inOTF->myDesc->jlocf]*
	     (inOTF->myDesc->crpix[inOTF->myDesc->jlocf]-1.0));
  nchan = inOTF->myDesc->inaxes[inOTF->myDesc->jlocf]; /* Number of freq channels */

  /* Average TRx, calJy */
  for (i=loDet; i<=hiDet; i++) {
    /* Channel flux density */
    ichan = i*inOTF->myDesc->incdatawt/inOTF->myDesc->incs;
    ichan = (ichan % nchan);
    ichan ++;
    freqRatio = ((inOTF->myDesc->crval[inOTF->myDesc->jlocf] + 
		 inOTF->myDesc->cdelt[inOTF->myDesc->jlocf]*
		 (ichan-inOTF->myDesc->crpix[inOTF->myDesc->jlocf])) /
		 refFreq);
    chFlux = FluxCal * pow(freqRatio,calIndex);
    if ((avgOff[i]!=fblank) && (avgOn[i]!=fblank)) {
      delta  = 0.5 *(avgOn[i] - avgOff[i]);
    } else {
      delta  = fblank;
    }
    if ((delta!=fblank) && (avgCal[i]!=fblank)) {
      TRx[i]     = avgOff[i] / avgCal[i];
      calJy[i]   = chFlux / (delta / avgCal[i]);
    } else {
      TRx[i]     = fblank;
      calJy[i]   = fblank;
    }
  } /* end loop over detectors */

  /* Save values on inOTF */
  dim[0] = ndetect; dim[1] = 1;
  ObitInfoListPut (inOTF->info, "tRx",     OBIT_float, dim,  (gpointer)TRx, err);
  ObitInfoListPut (inOTF->info, "calJy",   OBIT_float, dim,  (gpointer)calJy, err);

  /* Cleanup */
  if (farr)    g_free(farr);
  if (calJy)   g_free(calJy);
  if (TRx)     g_free(TRx);
  if (avgCal)  g_free(avgCal);
  if (avgOn)   g_free(avgOn);
  if (avgOff)  g_free(avgOff);
 
}  /* end ObitOTFCalUtilFitNod */

/**
 * Fits tipping scan
 * This routine presumes selected scan covers a large range in elevation.
 * Input values on inOTF
 * \li "Scan"     OBIT_Int (1,1,1) Scan number to process
 * \li "tCal"     OBIT_float (*,1,1) Noise cal value in K, per detector
 * \li "minEl"    OBIT_float (1,1,1) Minimum elevation allowed (deg)
 * \li "tSky"     OBIT_float (1,1,1) Physical temperature of the sky (K) [def 290]
 * \li "doPlot"   OBIT_Bool (1,1,1) If present and true, plot data and models [def FALSE]
 * Output values on inOTF
 * \li "Tau0"     OBIT_float (1,1,1) Zenith opacity in nepers
 *                i.e. the additional offset per airmass due to the atmosphere.
 *                per detector in units of the cal.
 * \li "aTemp"    OBIT_float (*,1,1) Effective atmospheric temperature in cal units.
 * \li "tRx"      OBIT_float (*,1,1) Receiver temperature per detector in units of the cal
 * \li "tipRMS"   OBIT_float (*,1,1) RMS residual (K) per detector
 * \param inOTF   Input OTF data. Maximum of MAXTIP samples used.
 * \param err     Error stack, returns if not empty.
 */
void ObitOTFCalUtilFitTip (ObitOTF *inOTF, ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat tau0, minEl, *ATemp=NULL, *TCal=NULL, *TRx=NULL, *RMS=NULL;
  olong i, ndetect, scan;
  ofloat TSky, TCalTatm, count, fblank = ObitMagicF();
  gchar *routine  ="ObitOTFCalUtilFitTip";

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));

  /* How many detectors? */
  ndetect = inOTF->geom->numberDetect;

  /* Create arrays */
  TCal    = g_malloc0(ndetect*sizeof(float));
  ATemp   = g_malloc0(ndetect*sizeof(float));
  TRx     = g_malloc0(ndetect*sizeof(float));
  RMS     = g_malloc0(ndetect*sizeof(float));

  /* Get control parameters */
  /* scan - no default */
  ObitInfoListGet (inOTF->info, "Scan", &type, dim,  (gpointer)&scan, err);
 
  /* tCal - no default */
  ObitInfoListGet (inOTF->info, "tCal", &type, dim,  (gpointer)TCal, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
 
  /* TSky - default aTemp */
  TSky = ATEMP;
  ObitInfoListGetTest (inOTF->info, "tSky", &type, dim,  (gpointer)&TSky);
 
  /* minEl - 0.0 */
  minEl = 0.0;
  ObitInfoListGetTest (inOTF->info, "minEl", &type, dim,  (gpointer)&minEl);
 
  /* Read, fit scan */
  FitTipScan (inOTF, scan, minEl, TRx, ATemp, RMS, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
  
  /* Compute tau0 from sky brightness (ATemp), 
     the assumed physical temperature (TSky),
     and the equivalent temperature of the cal (TCal) */
  /* Average the TCal * Tatm for the detectors  */
  TCalTatm = 0.0;
  count    = 0.0;
  for (i=0; i<ndetect; i++) {
    if (RMS[i] != fblank) {
      TCalTatm += TCal[i] * ATemp[i];
      count += 1.0;
    }
  }

  /* Average valid */
  if (count>0.0) TCalTatm /= count;
  else TCalTatm = TSky;

  tau0 = TCalTatm / TSky;
  /* debug
  fprintf (stderr,"tau0 = %f\n", tau0); */

  /* Convert RMSes to K */
  for (i=0; i<ndetect; i++) RMS[i] *= TCal[i];
  
  /* Save values on inOTF */
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (inOTF->info, "Tau0", OBIT_float, dim, (gpointer)&tau0, err);
  dim[0] = ndetect; dim[1] = 1;
  ObitInfoListPut (inOTF->info, "aTemp",   OBIT_float, dim,  (gpointer)ATemp, err);
  ObitInfoListPut (inOTF->info, "tRx",     OBIT_float, dim,  (gpointer)TRx, err);
  ObitInfoListPut (inOTF->info, "tipRMS",  OBIT_float, dim,  (gpointer)RMS, err);
  
  /* Cleanup */
  if (ATemp)   g_free(ATemp);
  if (TCal)    g_free(TCal);
  if (TRx)     g_free(TRx);
  if (RMS)     g_free(RMS);
  
}  /* end ObitOTFCalUtilFitTip */

/**
 * Adds flagging entry to associated flag table
 * Input values on inOTF
 * \li "flagVer"   OBIT_Int (1,1,1) Flagging table version, default = 1
 * \li "Chans"     OBIT_Int (2,1,1) First and highest channels to flag (1-rel), def, 0=>all
 * \li "timeRange" OBIT_float (2,1,1) Start and stop times to flag (days) def, 0s=all
 * \li "Target"    OBIT_string (?,1,1) Name of target, def, "Any" => all
 * \li "Stokes"    OBIT_string (?,1,1) Stokes to flag, def " " = flag all
 *                 "FFF"  where F is '1' to flag corresponding stokes, '0' not.
 *                 Stokes order 'R', 'L', 'RL' or 'X', 'Y', 'XY'
 * \li "Reason"    OBIT_string (?,1,1) reason string for flagging (max. 24 char).
 * \param inOTF   Input OTF data
 * \param err     Error stack, returns if not empty.
 * \return IO return code, OBIT_IO_OK = OK
 */
ObitIOCode ObitOTFCalUtilFlag (ObitOTF *inOTF, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  oint iarr[2];
  olong feed, flagVer, chans[2], TargID;
  ofloat timerange[2];
  gchar target[49], stokes[10], reason[49];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ObitTableOTFFlag    *FlagTable=NULL;
  ObitTableOTFFlagRow *FlagRow=NULL;
  ObitTableOTFTarget  *TargetTable=NULL;
  gchar *tname;
  olong ver, iRow;
  gchar *routine = "ObitOTFCalUtilFlag";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitOTFIsA(inOTF));

  /* Get parameters */
  flagVer = 1;
  ObitInfoListGetTest(inOTF->info, "flagVer", &type, dim, &flagVer);

  feed = 1;
  ObitInfoListGetTest(inOTF->info, "Feed", &type, dim, &feed);

  chans[0] = 1; chans[0] = 0;
  ObitInfoListGetTest(inOTF->info, "Chans", &type, dim, chans);

  timerange[0] = -1.0e20;  timerange[1] = 1.0e20;
  ObitInfoListGetTest(inOTF->info, "timeRange", &type, dim, timerange);
  /* default */
  if ((timerange[0]==0.0) && (timerange[1]==0.0)) {
    timerange[0] = -1.0e20;
    timerange[1] =  1.0e20;
  }

  g_snprintf (target, 48, "Any");
  ObitInfoListGetTest(inOTF->info, "Target", &type, dim, target);
  target[dim[0]] = 0;   /* terminate */

  g_snprintf (stokes, 9, " ");
  ObitInfoListGetTest(inOTF->info, "Stokes", &type, dim, stokes);
  stokes[dim[0]] = 0;   /* terminate */

  g_snprintf (reason, 48, " ");
  ObitInfoListGetTest(inOTF->info, "Reason", &type, dim, reason);
  reason[dim[0]] = 0;   /* terminate */

  /* Open/close input OTF to fully instantiate */
  retCode = ObitOTFOpen (inOTF, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, retCode);
  
  /* Close */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, retCode);

  /* Look up Target number if needed */
  if (strncmp ("Any", target, 3)) {
    /* Instantiate/Create Target Table */
    retCode = OBIT_IO_ReadErr;
    tname = g_strconcat ("OTFTarget table for: ",inOTF->name, NULL);
    ver = 1;
    TargetTable = newObitTableOTFTargetValue(tname, (ObitData*)inOTF, &ver, 
					     OBIT_IO_ReadWrite, err);
    dim[0] = strlen(target);
    dim[1] = 1;
    ObitTableOTFTargetLookup (TargetTable, dim, target, iarr, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, retCode);
    TargID = iarr[0];  /* Target ID */
    TargetTable = ObitTableOTFTargetUnref(TargetTable);
    g_free (tname);
  } else { /* flag all targets */
    TargID = 0;
  }

  /* Instantiate/Create output Flag Table */
  tname = g_strconcat ("OTFFlag table for: ", inOTF->name, NULL);
  ver = flagVer;
  FlagTable = newObitTableOTFFlagValue(tname, (ObitData*)inOTF, &ver, 
				       OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, retCode);
  g_free (tname);

  /* Open table */
  retCode = ObitTableOTFFlagOpen (FlagTable, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, retCode);

  /* Create Table Row */
  FlagRow = newObitTableOTFFlagRow (FlagTable);
  
  /* Attach  row to output buffer */
  ObitTableOTFFlagSetRow (FlagTable, FlagRow, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, retCode);

  /* If there are entries in the table, mark it unsorted */
  if (FlagTable->myDesc->nrow>0) 
    {FlagTable->myDesc->sort[0]=0; FlagTable->myDesc->sort[1]=0;}
  
  /* Fill in Flag row */
  FlagRow->TargetID = TargID;
  FlagRow->Feed     = feed;
  FlagRow->TimeRange[0] = timerange[0];
  FlagRow->TimeRange[1] = timerange[1];
  FlagRow->chans[0] = chans[0];
  FlagRow->chans[1] = chans[1];
  FlagRow->pFlags[0] = 0;
  strncpy (FlagRow->reason, reason, 24);
  if (stokes[0]==' ') {
    FlagRow->pFlags[0] = 15;
  } else {
    if (stokes[0]!='0') FlagRow->pFlags[0] += 1;
    if (stokes[1]!='0') FlagRow->pFlags[0] += 2;
    if (stokes[2]!='0') FlagRow->pFlags[0] += 4;
  }
  
  /* write row */
  iRow = FlagTable->myDesc->nrow+1;
  retCode =  ObitTableOTFFlagWriteRow (FlagTable, iRow, FlagRow, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, retCode);
   
  /* Close Flag table */
  retCode =  ObitTableOTFFlagClose (FlagTable, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, retCode);

  /* Cleanup */
  FlagTable = ObitTableOTFFlagUnref(FlagTable);
  FlagRow   = ObitTableOTFFlagRowUnref(FlagRow);

  return retCode;
} /* end ObitOTFCalUtilFlag */

/*---------------Private functions--------------------------*/

/**
 * Fits single calibrator scan.
 * This routines presumes that the selected detector(s) is scanned linearly through the 
 * expected position of the calibrator .
 * Gets calibrator information from the Target table (position, flux density)
 * \param inOTF    Input OTF data. 
 * \param detect   detector number (0-rel), -1 => all
 * \param scan     which scan to process
 * \param TRx      [out] Receiver (+atm) sky value in units of the cal per detector
 *                 set to fblank if no data.
 * \param calJy    [out] The calibrator in Jy, per detector
 * \param RAoff    [out] The RA offset (deg) per detector
 * \param Decoff   [out] The Dec offset (deg) per detector
 * \param Timeoff  [out] The Time offset (day) per detector
 * \param data     Work array MAXSAMPLE*ndetect in size
 * \param elev     [out] elevation (deg) at center of scan 
 * \param err      Error stack, returns if not empty.
 */
 static void FitCalScan (ObitOTF *inOTF, olong detect, olong scan,
			 ofloat *TRx, ofloat *calJy, 
			 ofloat *RAoff, ofloat *Decoff, ofloat *Timeoff,
			 ofloat *data, ofloat *elev, ObitErr *err)
  {
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOCode retCode;
  ObitInfoType type;
  ObitTableOTFTarget* targetTable=NULL;
  gboolean doCalSelect, doPlot, someOK, someData=FALSE;
  olong i, doCal, scans[2], ndetect,  mdetect, loDet, hiDet, iDet, targID=0;
  olong incdatawt;
  olong Count;
  ofloat Time[MAXSAMPLE], RA[MAXSAMPLE], Dec[MAXSAMPLE], work[MAXSAMPLE];
  gboolean isCal[MAXSAMPLE];
  ofloat *rec, fblank = ObitMagicF();
  ofloat peak, val, base, cal, center, sigma, FWHM;
  odouble sum, sum1, sum2, sumBase, mom1, mom2;
  odouble sumRA, sumDec, sumTime;
  odouble RACal, DecCal;
  ofloat FluxCal, ppeak, step;
  olong cntBase, cnt, peakpos=0, ver;
  gchar *routine = "FitCalScan";

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));

  /* How many detectors? */
  ndetect = inOTF->geom->numberDetect;
  /* How many to use */
  if (detect<0) {
    mdetect = ndetect;
    loDet = 0;
    hiDet = ndetect-1;
  } else {
    mdetect = 1;
    loDet = detect;
    hiDet = detect;
  }

  /* Want plotting? */
  doPlot = FALSE;
  ObitInfoListGetTest (inOTF->info, "doPlot", &type, dim,  (gpointer)&doPlot);

 /* Select scan on input */
  doCalSelect = TRUE;
  dim[0] = 1;
  ObitInfoListAlwaysPut(inOTF->info, "doCalSelect", OBIT_bool, dim, &doCalSelect);
  doCal = 1;
  ObitInfoListAlwaysPut(inOTF->info, "doCalib", OBIT_bool, dim, &doCal);
  scans[0] = scan; scans[1] = scan;
  dim[0] = 2;
  ObitInfoListAlwaysPut(inOTF->info, "Scans", OBIT_long, dim, scans);
  incdatawt = inOTF->myDesc->incdatawt; /* increment in data-wt axis */

   /* open OTF data to fully instantiate  */
  retCode = ObitOTFOpen (inOTF, OBIT_IO_ReadCal, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Read data into storage */
  Count = 0;      /* no. values in arrays */

  /* loop reading data */
  retCode = OBIT_IO_OK;
  while (retCode == OBIT_IO_OK) {

    /* read buffer */
    retCode = ObitOTFReadSelect (inOTF, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
    if (retCode==OBIT_IO_EOF) break; /* done? */

    /* Record pointer */
    rec = inOTF->buffer;
  
    /* Loop over buffer */
       for (i=0; i<inOTF->myDesc->numRecBuff; i++) {

	 /* Check -debug
	 if (scan!=((olong)(rec[inOTF->myDesc->ilocscan]+0.5))) {
	   fprintf (stderr,"Holy shit Batman this %d is not scan %d\n",
		     ((olong)(rec[inOTF->myDesc->ilocscan]+0.5)), scan);
	   exit(9);
	 } */
	   
	 /* Accumulate values  */
	 if (Count<MAXSAMPLE) {	   
	   isCal[Count] = rec[inOTF->myDesc->iloccal]!=0.0; /* Cal on? */
	   Time[Count]  = rec[inOTF->myDesc->iloct];
	   RA[Count]    = rec[inOTF->myDesc->ilocra];
	   Dec[Count]   = rec[inOTF->myDesc->ilocdec];
	   targID       = (olong)rec[inOTF->myDesc->iloctar];
	   for (iDet=loDet; iDet<=hiDet; iDet++) 
	     data[(iDet-loDet)*MAXSAMPLE+Count] = 
	       rec[inOTF->myDesc->ilocdata+iDet*incdatawt];
	   Count++;
	 } /* end  accumulate */
	 rec += inOTF->myDesc->lrec; /* Data record pointer */	 
       } /* end loop over buffer load */
  } /* end loop reading data */ 
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Better have some data */
  if (Count<=0) {
    Obit_log_error(err, OBIT_Error, "%s: NO data selected in %s", 
		   routine, inOTF->name);
    return;
 }

  /* Get elevation at center */
  *elev = ObitOTFArrayGeomElev (inOTF->geom, Time[Count/2], RA[Count/2], Dec[Count/2]);

   /* Get source info */
  ver = 1;
  targetTable = 
    newObitTableOTFTargetValue ("TargetTable", (ObitData*)inOTF, &ver, OBIT_IO_ReadWrite, 
				err);
  ObitTableOTFTargetGetSource (targetTable, targID, &RACal, &DecCal, &FluxCal, err);
  targetTable = ObitTableOTFTargetUnref(targetTable);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Make sure there is something */
  if ((RACal==0.0) || (DecCal==0.0) || (FluxCal==0.0)) {
    Obit_log_error(err, OBIT_Error, "%s: MISSING Calibrator info for %d %lf %lf %f in %s", 
		   routine, targID, RACal, DecCal, FluxCal, inOTF->name);
  }

  /* convert time to seconds */
  for (i=0;  i<Count; i++) Time[i] *= 86400.0;
      
  /* Convert positions and time to those relative to expected passage 
     through the source position */
  SetWrtCal (RACal, DecCal, Count, RA, Dec, Time);

  /* Average angle step per time unit */
  step = sqrt ((RA[Count-1]-RA[0])*(RA[Count-1]-RA[0])*cos(DecCal*DG2RAD)*cos(DecCal*DG2RAD) + 
	       (Dec[Count-1]-Dec[0])*(Dec[Count-1]-Dec[0])) / (Time[Count-1]-Time[0]);

  /* Loop over detectors */
  for (iDet=0; iDet<mdetect; iDet++) {

    /* Copy to work array - ObitOTFGetSolnAvgCal clobbers input */
    for (i=0;  i<Count; i++) work[i] = data[iDet*MAXSAMPLE+i];
    /* get average cal value */
    cal = ObitOTFGetSolnAvgCal (Count, isCal, work);
      
    /* Get basic statistics, moments */
    someOK   = FALSE;
    calJy[iDet+loDet] = fblank;  /* in case no data */
    TRx[iDet+loDet]   = fblank;  /* in case no data */
    peak = -1.0E20;
    sumBase = 0.0;
    cntBase = 0;
    for (i=0; i<Count; i++) {

      /* valid? */
      val = data[iDet*MAXSAMPLE+i];
      if (data[iDet*MAXSAMPLE+i]==fblank) continue;
      someOK   = TRUE;
      someData = TRUE;

      /* If cal on, correct */
      if (isCal[i]) data[iDet*MAXSAMPLE+i] -= cal;
      val = data[iDet*MAXSAMPLE+i];

      /* Use first and last 10% as baseline */
      if ((i<(0.1*Count)) || (i>(0.9*Count))) {
	sumBase += val;
	cntBase++;
      }

      /* Look for peak */
      if (val>peak) {
	peak = val;
	peakpos = i;
      }

    } /* End first loop over data */
    ppeak = peak;  /* save for diagnostics */

    /* Find any good data? */
    if (!someOK) continue;

    /* Baseline value */
    if (cntBase>0) base = sumBase / cntBase;
    else base = 0.0;
    peak -= base;  /* want height of peak above baseline */

    /* Get first moment */
    sum = sum1 = 0.0;
    sumRA = sumDec = sumTime = 0.0;
    cnt = 0;
    for (i=0; i<Count; i++) {
      /* valid? */
      val = data[iDet*MAXSAMPLE+i];
      if (val==fblank) continue; /* ignore blanked data */
      val -= base;    /* correct for baseline*/
      sum  += val;
      sum1 += ((odouble)i) * val;
      sumRA   += val * RA[i];
      sumDec  += val * Dec[i];
      sumTime += val * Time[i];
      cnt++;
    } /* End loop over data for first moment */
    /* Get first moments */
    if (sum>0) {
      RAoff[iDet+loDet]   = (ofloat)(sumRA  / sum);  /* RA offset */
      Decoff[iDet+loDet]  = (ofloat)(sumDec  / sum); /* Dec offset */
      Timeoff[iDet+loDet] = (ofloat)(sumTime / sum); /* Time offset */
      mom1 = sum1/sum;
    } else {  /* Use position of peak */
      RAoff[iDet+loDet]   = RA[peakpos];
      Decoff[iDet+loDet]  = Dec[peakpos];
      Timeoff[iDet+loDet] = Time[peakpos];
      mom1 = peakpos;
    }

    /* Get second moment */
    sum = sum2 =0.0;
    cnt = 0;
    for (i=0; i<Count; i++) {
      /* valid? */
      val = data[iDet*MAXSAMPLE+i];
      if (val==fblank) continue;
      val -= base;   /* correct for baseline*/
      sum += val;
      sum2 += ((odouble)i - mom1) * ((odouble)i - mom1) * val;
      cnt++;
    } /* End loop over data for second moment */
    if (sum>0.0)  mom2 = sqrt (fabs(sum2) / sum);
    else mom2 = 2.0;

    /* Do ls fit here to Gaussian + baseline */
    center = Timeoff[iDet+loDet];
    sigma = 15.0 * 0.1*fabs((Time[Count/2-5]-Time[Count/2+5]));
    FitCalGauss (Count, Time, &data[iDet*MAXSAMPLE], &base, &peak, &center, &sigma, err);
    if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

    /* sigma to degrees */
    sigma *= step;
    FWHM = sigma * 2.355 *3600.0; /* FWHM in asec */

    /* Plot if requested */
    if (doPlot)
      PlotCalGauss (Count, Time, scan, iDet+1, &data[iDet*MAXSAMPLE], 
		    base, peak, center, sigma, err);

    /* Get results */
    calJy[iDet+loDet] = FluxCal * cal / peak;
    TRx[iDet+loDet]   = base / cal;

    /* DEBUG */
    fprintf (stdout, "detector %d base %f peak %f pos %d cal %f flux %f RACal %lf DecCal %lf\n",
	    iDet+loDet, base, ppeak, peakpos, cal, FluxCal, RACal, DecCal );
    fprintf (stdout, "scan %d calJy %f TRx %f peak %f mom1 %f mom2 %f\n",
	     scan, calJy[iDet+loDet], TRx[iDet+loDet], peak, mom1, mom2);
    fprintf (stdout, "RAoff %lf Decoff %f Timeoff %f sigma %f FWHM %f\n", 
	     RAoff[iDet+loDet]*3600.0,Decoff[iDet+loDet]*3600.0, 
	     Timeoff[iDet+loDet],sigma, FWHM );
    

    /* Check for very bad position offsets (error) */
    if ((fabs(RAoff[iDet+loDet])>1.0) || (fabs(Decoff[iDet+loDet])>1.0)) {
      Obit_log_error(err, OBIT_InfoWarn, 
		     "Warning: Excessive offset in RA (%f) or Dec(%f) detector %d scan %d",
		     RAoff[iDet+loDet], Decoff[iDet+loDet], iDet, scan);
    }
    
  } /* end loop over detectors */

  /* Give warning if no data selected */
  if (!someData) Obit_log_error(err, OBIT_InfoWarn, 
				"%s: Warning: NO data selected", routine);

}  /* end FitCalScan */

/**
 * Averages the value in a given scan
 * \param inOTF    Input OTF data. 
 * \param detect   detector number (0-rel), -1 => all
 * \param scan     which scan to process
 * \param calOff   [out] Average value cal off
 * \param cal      [out] Average cal value
 * \param data     Work array MAXSAMPLE*ndetect in size,
 *                 Returned with scan data with Cal corrected.
 * \param targID   [out] Target ID of scan
 * \param err      Error stack, returns if not empty.
 */
 static void FitCalAverage (ObitOTF *inOTF, olong detect, olong scan,
			 ofloat *calOff, ofloat *cal, 
			 ofloat *data, olong *targID, ObitErr *err)
  {
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOCode retCode;
  gboolean doCalSelect;
  olong i, doCal, scans[2], ndetect,  mdetect, loDet, hiDet, iDet;
  olong incdatawt;
  olong Count;
  ofloat val, work[MAXSAMPLE];
  gboolean isCal[MAXSAMPLE];
  ofloat *rec, calAvg, fblank = ObitMagicF();
  odouble sum;
  olong cnt;
  gchar *routine = "FitCalAverage";

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));

  /* How many detectors? */
  ndetect = inOTF->geom->numberDetect;
  /* How many to use */
  if (detect<0) {
    mdetect = ndetect;
    loDet = 0;
    hiDet = ndetect-1;
  } else {
    mdetect = 1;
    loDet = detect;
    hiDet = detect;
  }

  /* Select scan on input */
  doCalSelect = TRUE;
  dim[0] = 1;
  ObitInfoListAlwaysPut(inOTF->info, "doCalSelect", OBIT_bool, dim, &doCalSelect);
  doCal = 1;
  ObitInfoListAlwaysPut(inOTF->info, "doCalib", OBIT_bool, dim, &doCal);
  scans[0] = scan; scans[1] = scan;
  dim[0] = 2;
  ObitInfoListAlwaysPut(inOTF->info, "Scans", OBIT_long, dim, scans);
  incdatawt = inOTF->myDesc->incdatawt; /* increment in data-wt axis */

   /* open OTF data to fully instantiate  */
  retCode = ObitOTFOpen (inOTF, OBIT_IO_ReadCal, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Read data into storage */
  Count = 0;      /* no. values in arrays */

  /* loop reading data */
  retCode = OBIT_IO_OK;
  while (retCode == OBIT_IO_OK) {

    /* read buffer */
    retCode = ObitOTFReadSelect (inOTF, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
    if (retCode==OBIT_IO_EOF) break; /* done? */

    /* Record pointer */
    rec = inOTF->buffer;
  
    /* Loop over buffer */
       for (i=0; i<inOTF->myDesc->numRecBuff; i++) {

	 /* Check -debug
	 if (scan!=((olong)(rec[inOTF->myDesc->ilocscan]+0.5))) {
	   fprintf (stderr,"Holy shit Batman this %d is not scan %d\n",
		     ((olong)(rec[inOTF->myDesc->ilocscan]+0.5)), scan);
	   exit(9);
	 } */
	   
	 /* Accumulate values  */
	 if (Count<MAXSAMPLE) {	   
	   isCal[Count] = rec[inOTF->myDesc->iloccal]!=0.0; /* Cal on? */
	   *targID      = (olong)rec[inOTF->myDesc->iloctar];
	   for (iDet=loDet; iDet<=hiDet; iDet++) 
	     data[(iDet-loDet)*MAXSAMPLE+Count] = 
	       rec[inOTF->myDesc->ilocdata+iDet*incdatawt];
	   Count++;
	 } /* end  accumulate */
	 rec += inOTF->myDesc->lrec; /* Data record pointer */	 
       } /* end loop over buffer load */
  } /* end loop reading data */ 
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Better have some data */
  if (Count<=0) {
    Obit_log_error(err, OBIT_Error, "%s: NO data selected in %s", 
		   routine, inOTF->name);
    return;
 }

  /* Loop over detectors averaging */
  for (iDet=0; iDet<mdetect; iDet++) {

    /* Copy to work array - ObitOTFGetSolnAvgCal clobbers input */
    for (i=0;  i<Count; i++) work[i] = data[iDet*MAXSAMPLE+i];
    /* get average cal value */
    calAvg = ObitOTFGetSolnAvgCal (Count, isCal, work);
      
    /* Average */
    sum = 0.0;
    cnt = 0;
    for (i=0; i<Count; i++) {

      /* valid? */
      val = data[iDet*MAXSAMPLE+i];
      if (val==fblank) continue;

      /* If cal on, correct */
      if (isCal[i]) data[iDet*MAXSAMPLE+i] -= calAvg;
      val = data[iDet*MAXSAMPLE+i];

      /* Sums */
      sum += val;
      cnt++;
    } /* End loop over data */

    /* save results */
    cal[iDet] = calAvg;
    if (cnt>0) calOff[iDet] = sum / cnt;
    else calOff[iDet] = fblank;
    /* debug 
    fprintf (stdout, "detector %d scan %d sum %f cnt %d calOff %f calAvg %f\n",
	    iDet+loDet, scan, sum, cnt, calOff[iDet], cal[iDet]);
   */
  } /* end loop over detectors */
}  /* end FitCalAverage */

/**
 * Averages the spectra in a given scan
 * \param inOTF    Input OTF data. 
 * \param scan     which scan to process
 * \param data     [out] Average selected spectra, fblanked
 * \param wt       [out] weight array the size of data
 * \param targID   [out] Target ID of scan
 * \param parAng   [out] Average parallactic angle (deg)
 * \param Time     [out] Average time  (day)
 * \param err      Error stack, returns if not empty.
 */
static void SpectrumAverage (ObitOTF *inOTF,olong scan, ofloat *data, 
			     ofloat *wt, olong *targID, ofloat *parAng, 
			     ofloat *Time, ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOCode retCode;
  gboolean doCalSelect;
  ObitOTFDesc *desc = inOTF->myDesc;
  olong i, doCal, scans[2], nfeed, nstok, nchan, Count, ifeed, istok, ichan, ndetect;
  olong incdatawt, incf, incs, incfeed, indx, jndx;
  ofloat *rec, sumPA, cntPA, sumTim, cntTim, fblank = ObitMagicF();
  gchar *routine = "SpectrumAverage";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));
  
  /* How many detectors? */
  nfeed   = desc->inaxes[desc->jlocfeed];
  nstok   = desc->inaxes[desc->jlocs];
   /* Stokes limited to 2 solutions */
  if (nstok>2) nstok = 2;
  nchan   = desc->inaxes[desc->jlocf];
  ndetect = nfeed * nstok * nchan;

  /* Select scan on input */
  doCalSelect = TRUE;
  dim[0] = 1;
  ObitInfoListAlwaysPut(inOTF->info, "doCalSelect", OBIT_bool, dim, &doCalSelect);
  doCal = 1;
  ObitInfoListAlwaysPut(inOTF->info, "doCalib", OBIT_bool, dim, &doCal);
  scans[0] = scan; scans[1] = scan;
  dim[0] = 2;
  ObitInfoListAlwaysPut(inOTF->info, "Scans", OBIT_long, dim, scans);
  incdatawt = desc->incdatawt; /* increment in data-wt axis */
  incf      = desc->incf;
  incs      = desc->incs;
  incfeed   = desc->incfeed;
 
  /* open OTF data  */
  retCode = ObitOTFOpen (inOTF, OBIT_IO_ReadCal, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
  
  /* Read data into storage */
  /* Init */
  for (i=0; i<ndetect; i++) data[i] = wt[i] = 0.0;
  Count = 0;  
  sumPA  = cntPA  = 0.0;
  sumTim = cntTim = 0.0;
  /* loop accumulating data */
  retCode = OBIT_IO_OK;
  while (retCode == OBIT_IO_OK) {
    
    /* read buffer */
    retCode = ObitOTFReadSelect (inOTF, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
    if (retCode==OBIT_IO_EOF) break; /* done? */
    
    /* Record pointer */
    rec = inOTF->buffer;
    
    /* Loop over buffer */
    for (i=0; i<inOTF->myDesc->numRecBuff; i++) {
      
      /* Accumulate values, order chan/Stokes/feed  */
      *targID  = (olong)rec[inOTF->myDesc->iloctar];
      sumPA    += rec[inOTF->myDesc->ilocrot];
      cntPA++;
      sumTim   += rec[inOTF->myDesc->iloct];
      cntTim++;
      jndx = 0;
      /* Loop over feed */
      for (ifeed=0; ifeed<nfeed; ifeed++) {
	/* Loop over Stokes */
	for (istok=0; istok<nstok; istok++) {
	  /* Loop over Channel */
	  for (ichan=0; ichan<nchan; ichan++) {
	    indx = desc->ilocdata + ichan*incf + istok*incs + ifeed*incfeed;
	    if ((rec[indx+1]>0.0) && (rec[indx]!=fblank)) {
	      data[jndx] += rec[indx]*rec[indx+1];
	      wt[jndx]   += rec[indx+1];
	    }
	    jndx++;
	    Count++;
	  } /* end channel loop */
	} /* end Stokes loop */
      } /* end feed loop */
      rec += inOTF->myDesc->lrec; /* Data record pointer */	 
    } /* end loop over buffer load */
  } /* end loop reading data */ 
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
  
  /* Better have some data */
  if (Count<=0) {
    Obit_log_error(err, OBIT_Error, "%s: NO data selected in %s", 
		   routine, inOTF->name);
    return;
  }
  
  /* Loop over detectors averaging */
  for (i=0; i<ndetect; i++) {
    if (wt[i]>0.0) data[i] /= wt[i];
    else           data[i] = fblank;
  } /* end averaging loop */

  /* Average parallactic angle */
  *parAng = sumPA/(cntPA+1.0e-9);
  *Time   = sumTim/(cntTim+1.0e-9);
}  /* end SpectrumAverage */

/**
 * Averages on source, off source and the and cal value for dual beam
 * switched data.
 * \param inOTF    Input OTF data. 
 * \param detect   detector number (0-rel), -1 => all
 * \param scan     which scan to process
 * \param avgOff   [out] Average value off source
 * \param avgOn    [out] Average value on source
 * \param avgCal   [out] Average cal value
 * \param targId   [out] Target ID of scan
 * \param err      Error stack, returns if not empty.
 */
static void FitCalNod (ObitOTF *inOTF, olong detect, olong scan,
		       ofloat *avgOff, ofloat *avgOn, ofloat *avgCal, 
		       olong *targID, ObitErr *err)
  {
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOCode retCode;
  ObitTableOTFTarget* targetTable=NULL;
  gboolean doCalSelect;
  olong i, itemp, doCal, scans[2], incfeed, ndetect,  mdetect, loDet, hiDet, iDet;
  olong incdatawt;
  olong Count, ver;
  ofloat val, *feedRA=NULL, *feedDec=NULL;
  olong *state=NULL;
  gboolean gotTarinfo=FALSE, isRef, *isCal=NULL;
  odouble RACal, DecCal;
  ofloat *rec, *data=NULL, FluxCal, dra, ddec, fblank = ObitMagicF();
  odouble sumOn, sumOff;
  olong cntOn, cntOff, nrec=0;
  gchar *routine = "FitCalNod";

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));

  /* How many detectors? */
  ndetect = inOTF->geom->numberDetect;
  /* How many to use */
  if (detect<0) {
    mdetect = ndetect;
    loDet = 0;
    hiDet = ndetect-1;
  } else {
    mdetect = 1;
    loDet = detect;
    hiDet = detect;
  }

  /* Select scan on input */
  doCalSelect = TRUE;
  dim[0] = 1;
  ObitInfoListAlwaysPut(inOTF->info, "doCalSelect", OBIT_bool, dim, &doCalSelect);
  doCal = 1;
  ObitInfoListAlwaysPut(inOTF->info, "doCalib", OBIT_bool, dim, &doCal);
  scans[0] = scan; scans[1] = scan;
  dim[0] = 2;
  ObitInfoListAlwaysPut(inOTF->info, "Scans", OBIT_long, dim, scans);
  incdatawt = inOTF->myDesc->incdatawt; /* increment in data-wt axis */

   /* open OTF data to fully instantiate  */
  retCode = ObitOTFOpen (inOTF, OBIT_IO_ReadCal, err);
  if (err->error) goto cleanup;

  /* Read data into storage */
  Count = 0;      /* no. values in arrays */

  /* loop reading data */
  retCode = OBIT_IO_OK;
  while (retCode == OBIT_IO_OK) {

    /* read buffer */
    retCode = ObitOTFReadSelect (inOTF, NULL, err);
    if (err->error) goto cleanup;
    if (retCode==OBIT_IO_EOF) break; /* done? */

    /* Record pointer */
    rec = inOTF->buffer;
  
    /* Loop over buffer */
       for (i=0; i<inOTF->myDesc->numRecBuff; i++) {

	 /* Get source info, create buffers first record */
	 if (!gotTarinfo) {
	   gotTarinfo = TRUE;
	   /* Get position from table */
	   ver = 1;
	   *targID     = (olong)rec[inOTF->myDesc->iloctar];
	   targetTable = 
	     newObitTableOTFTargetValue ("TargetTable", (ObitData*)inOTF, &ver, OBIT_IO_ReadWrite, 
					 err);
	   ObitTableOTFTargetGetSource (targetTable, *targID, &RACal, &DecCal, &FluxCal, err);
	   targetTable = ObitTableOTFTargetUnref(targetTable);
	   if (err->error) goto cleanup;
	   
	   /* Make sure there is something */
	   if ((RACal==0.0) || (DecCal==0.0)) {
	     Obit_log_error(err, OBIT_Error, "%s: MISSING Calibrator info for %d %lf %lf in %s", 
			    routine, *targID, RACal, DecCal, inOTF->name);
	     goto cleanup;
	   }
	   /* Create arrays */
	   nrec = ObitOTFNumRecScan (inOTF);
	   feedRA  = g_malloc(ndetect*sizeof(ofloat));
	   feedDec = g_malloc(ndetect*sizeof(ofloat));
	   isCal   = g_malloc(nrec*sizeof(gboolean));
	   state   = g_malloc(nrec*sizeof(olong));
	   data    = g_malloc(ndetect*nrec*sizeof(ofloat));

	 } /* End of get target info & create arrays */

	 /* Get feed positions */
	 ObitOTFArrayGeomCoord (inOTF->geom,  rec[inOTF->myDesc->ilocra], 
				rec[inOTF->myDesc->ilocdec], rec[inOTF->myDesc->ilocrot], 
				feedRA, feedDec);
	 /* Three states here, sig beam on source (state=1), ref beam on source (state=-1) 
	    or neither (state=0) */
	 state[Count] = 0;
	 /* Close to sig beam (0)? */
	 itemp = 0;
	 dra  = feedRA[itemp]  - RACal;
	 ddec = feedDec[itemp] - DecCal;
	 if (sqrt(dra*dra+ddec*ddec) < 0.3*inOTF->myDesc->beamSize) state[Count] = 1;
	 if (!state[Count]) {
	   /* Close to reference beam? */
	   itemp = inOTF->myDesc->incfeed / inOTF->myDesc->incdatawt;
	   dra  = feedRA[itemp]  - RACal;
	   ddec = feedDec[itemp] - DecCal;
	   if (sqrt(dra*dra+ddec*ddec) < 0.3*inOTF->myDesc->beamSize) state[Count] = -1;
	   /* debug if (state[Count]==-1) {
	      fprintf (stderr,"Holy Shit Batman \n");
	      }*/
	 }

	 /* Data gives sig beam position, need ref beam */

	 /* Check -debug
	 if (scan!=((olong)(rec[inOTF->myDesc->ilocscan]+0.5))) {
	   fprintf (stderr,"Holy shit Batman this %d is not scan %d\n",
		     ((olong)(rec[inOTF->myDesc->ilocscan]+0.5)), scan);
	   exit(9);
	 } */
	   
	 /* Accumulate values  */
	 isCal[Count] = (state[Count]!=0) && rec[inOTF->myDesc->iloccal]!=0.0; /* Cal on? */
	 for (iDet=loDet; iDet<=hiDet; iDet++) 
	   data[(iDet-loDet)*nrec+Count] = 
	     rec[inOTF->myDesc->ilocdata+iDet*incdatawt];
	 Count++;
	 rec += inOTF->myDesc->lrec; /* Data record pointer */	 
       } /* end loop over buffer load */
  } /* end loop reading data */ 
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) goto cleanup;

  /* Better have some data */
  if (Count<=0) {
    Obit_log_error(err, OBIT_Error, "%s: NO data selected in %s", 
		   routine, inOTF->name);
    goto cleanup;
 }

  /* Loop over detectors averaging */
  for (iDet=0; iDet<mdetect; iDet++) {

    /* Is this a sig or ref beam */
    if (inOTF->myDesc->jlocfeed>=0) 
      incfeed = inOTF->myDesc->incfeed / inOTF->myDesc->incdatawt;
    else incfeed = 1;  /* This is probably a bad sign */
    itemp = iDet / incfeed;
    /* itemp odd is reference beam */
    isRef = itemp != 2*(itemp/2);

    /* Average */
    sumOn  = 0.0;
    cntOn  = 0;
    sumOff = 0.0;
    cntOff = 0;
    for (i=0; i<Count; i++) {

      /* valid? */
      val = data[iDet*nrec+i];
      if ((val==fblank) || (state[i]==0) || isCal[i]) continue;

      /* Sums */
      if (isRef) { /* reference beam */
	if (state[i]<0)      {sumOn  += val; cntOn++;}
	else if (state[i]>0) {sumOff += val; cntOff++;}
      } else { /* signal beam */
	if (state[i]>0)      {sumOn  += val; cntOn++;}
	else if (state[i]<0) {sumOff += val; cntOff++;}
      }
    } /* End loop over data */

    /* save results */
    if (cntOn>0)  avgOn[iDet]  = sumOn  / cntOn;
    if (cntOff>0) avgOff[iDet] = sumOff / cntOff;

    /* get average cal value */
    avgCal[iDet] = fabs (ObitOTFGetSolnAvgCal (Count, isCal, &data[iDet*nrec]));
      
    /* debug 
    fprintf (stdout, "detector %d scan %d sum %f cnt %d calOff %f calAvg %f\n",
	    iDet+loDet, scan, sum, cnt, calOff[iDet], cal[iDet]);
   */
  } /* end loop over detectors */

  cleanup:
  if (feedRA)  g_free(feedRA);
  if (feedDec) g_free(feedDec);
  if (isCal)   g_free(isCal);
  if (state)   g_free(state);
  if (data)    g_free(data);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
}  /* end FitCalNod */


/**
 * Given arrays of RA, Dec and time, all presumed to be sampled at the same times,
 * convert the position arrays to those relative to the source position and 
 * the corresponding time array to times relative expected passage through 
 * the source position.
 * \param RACal  Source RA
 * \param DecCal Source Dec
 * \param Count
 * \param RA 
 * \param Dec
 * \param Time
 */
static void SetWrtCal (odouble RACal, odouble DecCal,
		       olong Count, ofloat *RA, ofloat *Dec, ofloat *Time)
{
  olong i,lo, hi, cnt;
  gboolean isRA;
  ofloat *pos;
  odouble sumX, sumXX, sumY, sumXY, a, b, T0;

  /* Convert positions to relative */
  for (i=0; i<Count; i++) {
    RA[i]  -= RACal;
    Dec[i] -= DecCal;
  }

  /* presume scan in either RA or Dec - decide which */
  isRA = (fabs(RA[Count-1]-RA[0]) > fabs(Dec[Count-1]-Dec[0]));
  if (isRA) pos = RA; /* Get pointer for direction being scanned */
  else pos = Dec;

  /* Determine zero from inner half of scan */
  lo = (olong)(0.25*Count+0.5);
  hi = (olong)(0.75*Count+0.5);
  sumX = sumXX = sumY = sumXY = 0.0;
  cnt = 0;
  for (i=lo; i<=hi; i++) {
    cnt++;
    sumY  += pos[i];
    sumX  += Time[i];
    sumXY += pos[i] * Time[i];
    sumXX += Time[i] * Time[i];
  }

  /* Linear regression coefficients */
  a = (sumY*sumXX - sumX*sumXY) / (cnt*sumXX - sumX*sumX);
  b = (cnt*sumXY  - sumX*sumY)  / (cnt*sumXX - sumX*sumX);

  /* Time of source passage */
  T0 = -a / b;

  /* Convert time to relative to source passage */
  for (i=0; i<Count; i++) Time[i] -= T0;
  
} /* end SetWrtCal */


/* Structure for least squares fitting */
struct fitData {
  size_t n;   /* number of data points */
  float *x;   /* independent variable */
  float *y;   /* dependent variable */
};

/**
 * Gaussian model calculating routine for gsl least squares fitter
 * \param coef   Coefficient array (baseline, amplitude, center, 1/variance)
 * \param params Data structure
 * \param f      [out] function residuals
 * \returns GSL completion code
 */
static int gaussFunc (const gsl_vector *coef, void *params, gsl_vector *f)
{
  size_t n = ((struct fitData *)params)->n;
  float *x = ((struct fitData *)params)->x;
  float *y = ((struct fitData *)params)->y;
  long i;
  double base, amp, cen, ivar, model, resid;

  /* Current parameters */
  base = gsl_vector_get(coef, 0);
  amp  = gsl_vector_get(coef, 1);
  cen  = gsl_vector_get(coef, 2);
  ivar  = gsl_vector_get(coef, 3);

  /* Loop through data calculating residuals to model */
  for (i=0; i<n; i++) {
    model = base + amp * exp (-(x[i]-cen)*(x[i]-cen) * ivar);
    resid = model - y[i];
    gsl_vector_set(f, i, resid);  /* to output vector */
  }

  return GSL_SUCCESS;
} /* end  gaussFunc */

/**
 * Gaussian model calculating Jacobean for gsl least squares fitter
 * This is the partial derivative of the residuals matrix
 * \param coef   Coefficient array (baseline, amplitude, center, 1/variance)
 * \param params Data structure
 * \param J      [out] Jacobean values
 * \returns GSL completion code
 */
static int gaussJacob (const gsl_vector *coef, void *params, gsl_matrix *J)
{
  size_t n = ((struct fitData *)params)->n;
  float *x = ((struct fitData *)params)->x;
  long i;
  double base, amp, cen, ivar, eterm;
  double part1, part2, part3, part4;

  /* Current parameters */
  base = gsl_vector_get(coef, 0);
  amp  = gsl_vector_get(coef, 1);
  cen  = gsl_vector_get(coef, 2);
  ivar = gsl_vector_get(coef, 3);

  /* Loop through data calculating partial derivatives of residuals */
  for (i=0; i<n; i++) {
    eterm = exp (-(x[i]-cen)*(x[i]-cen) * ivar);
    part1 = 1.0;                              /* partial wrt baseline */
    part2 = eterm;                            /* partial wrt amplitude */
    part3 = +amp*eterm*ivar*2.0*(x[i]-cen);   /* partial wrt center */
    part4 = -amp*eterm*(x[i]-cen)*(x[i]-cen); /* partial wrt 1/var  */
    gsl_matrix_set(J, i, 0, part1);  /* to output matrix */
    gsl_matrix_set(J, i, 1, part2);
    gsl_matrix_set(J, i, 2, part3);
    gsl_matrix_set(J, i, 3, part4);
  }

  return GSL_SUCCESS;
} /* end  gaussJacob */


/**
 * Compute both function and derivatives for  gsl least squares fitter
 * \param coef   Coefficient array (baseline, amplitude, center, sigma)
 * \param params Data structure
 * \param f      [out] function residuals
 * \param J      [out] Jacobean
 * \returns GSL completion code
 */
static int gaussFuncJacob (const gsl_vector *coef, void *params, gsl_vector *f, 
			    gsl_matrix *J)
{
  gaussFunc(coef, params, f);
  gaussJacob(coef, params, J);
  return GSL_SUCCESS;
} /* end  gaussFuncJacob */

/**
 * Fit Gaussian + baseline to a time sequence of data
 * \param Count  Number of samples to be fitted
 * \param Time   Time tags of the data 
 * \param Temp   Measured values
 * \param base   [in/out] value off source
 * \param peak   [in/out] peak value
 * \param center [in/out] Time of center of Gaussian
 * \param sigma  [in/out] Gaussian sigma (width) in time units.
 * \param err      Error stack, returns if not empty.
 */
static void FitCalGauss (olong Count, ofloat *Time, ofloat *Temp,  
			 ofloat *base, ofloat *peak, ofloat *center, 
			 ofloat *sigma, ObitErr *err)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver* solver;
  gsl_multifit_function_fdf func;
  gsl_vector_view coef;
  struct fitData data;
  int status, iter;
  double coef_init[4];
  gchar *routine = "FitCalGauss";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  if (Count<=0) return;

  /* Fitter data structure */
  data.n = Count;
  data.x = Time;
  data.y = Temp;

  /* initial guess */
  coef_init[0] = *base;
  coef_init[1] = *peak;
  coef_init[2] = *center;
  coef_init[3] = 1.0 / (*sigma * *sigma); /* fit as 1/var */

  coef = gsl_vector_view_array (coef_init, 4);

  /* Create /fill function structure */
  func.f   = &gaussFunc;      /* Compute function */
  func.df  = &gaussJacob;     /* Compute Jacobian (derivative matrix) */
  func.fdf = &gaussFuncJacob; /* Compute both function and derivatives */
  func.n = Count;             /* Number of data points */
  func.p = 4;                 /* number of parameters */
  func.params = &data;        /* Data structure */

  T = gsl_multifit_fdfsolver_lmsder;
  solver = gsl_multifit_fdfsolver_alloc(T, Count, 4);
  gsl_multifit_fdfsolver_set(solver, &func, &coef.vector);

  /* ready to rumble */
  iter = 0;
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);
    if ((status!=GSL_CONTINUE) && (status!=GSL_SUCCESS)) {/* problem? */
      Obit_log_error(err, OBIT_Error, "%s: Solver status %d %s", 
		     routine,status,gsl_strerror(status));
      break; 
    }
    /* convergence test */
    status = gsl_multifit_test_delta(solver->dx, solver->x, 1.0e-4, 1.0e-4);
   }
  while ((status == GSL_CONTINUE) && (iter< 500));

  /* return results */
  /* debug fprintf (stderr,"no. iter %d base %f peak %f\n",iter, *base,*peak);*/
  
  *base   = (float)gsl_vector_get(solver->x, 0);
  *peak   = (float)gsl_vector_get(solver->x, 1);
  *center = (float)gsl_vector_get(solver->x, 2);
  *sigma  = (float)gsl_vector_get(solver->x, 3);
  *sigma = 1.0 / (sqrt(*sigma)); /* back to sigma */

  /* cleanup */
  gsl_multifit_fdfsolver_free(solver);
} /* end FitCalGauss */


/**
 * Fits tipping scan to determine reciever temperature and sky brightness
 * This routine presumes selected scan covers a large range in elevation.
 * Assumes sky brightness increases linearly with air mass (breaks down 
 * for large opacity)
 * \param inOTF    Input OTF data. 
 * \param scan     which scan to process
 * \param minEl    minimum elevation (deg) to use.
 * \param TRx      [out] Reciever (+atm) sky value in units of the cal per detector
 * \param ATemp    [out] The sky brightness in cal units per airmass
 * \param RMS      [out] The RMS residual in cal units, fblank if undetermined
 * \param err      Error stack, returns if not empty.
 */
static void FitTipScan (ObitOTF *inOTF, olong scan,  ofloat minEl,
			ofloat *TRx, ofloat *ATemp, ofloat *RMS, ObitErr *err)
  {
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOCode retCode;
  ObitInfoType type;
  gboolean doCalSelect, doPlot;
  olong i, j, doCal, scans[2], ndetect, incdatawt;
  olong Count;
  ofloat airmass[MAXTIP], elev, work[MAXTIP], *data=NULL;
  gboolean isCal[MAXTIP];
  ofloat *rec, *TCal=NULL, minElR, fblank = ObitMagicF();
  ofloat TSky, tau0, cal, ical;
  gchar *routine = "FitTipScan";

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));

  /* How many detectors? */
  ndetect = inOTF->geom->numberDetect;

  /* allocate data array */
  data = g_malloc0(MAXTIP*ndetect*sizeof(float));
  minElR = minEl*1.745329e-2; /* to radians */

  /* Select scan on input */
  doCalSelect = TRUE;
  dim[0] = 1;
  ObitInfoListAlwaysPut(inOTF->info, "doCalSelect", OBIT_bool, dim, &doCalSelect);
  doCal = 1;
  ObitInfoListAlwaysPut(inOTF->info, "doCalib", OBIT_bool, dim, &doCal);
  scans[0] = scan; scans[1] = scan;
  dim[0] = 2;
  ObitInfoListAlwaysPut(inOTF->info, "Scans", OBIT_long, dim, scans);
  incdatawt = inOTF->myDesc->incdatawt; /* increment in data-wt axis */

   /* open OTF data to fully instantiate  */
  retCode = ObitOTFOpen (inOTF, OBIT_IO_ReadCal, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Read data into storage */
  Count = 0;      /* no. values in arrays */

  /* loop reading data */
  retCode = OBIT_IO_OK;
  while (retCode == OBIT_IO_OK) {

    /* read buffer */
    retCode = ObitOTFReadSelect (inOTF, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
    if (retCode==OBIT_IO_EOF) break; /* done? */

    /* Record pointer */
    rec = inOTF->buffer;
  
    /* Loop over buffer */
       for (i=0; i<inOTF->myDesc->numRecBuff; i++) {

	 /* Check -debug
	 if (scan!=((olong)(rec[inOTF->myDesc->ilocscan]+0.5))) {
	   fprintf (stderr,"Holy shit Batman this %d is not scan %d\n",
		     ((olong)(rec[inOTF->myDesc->ilocscan]+0.5)), scan);
	   exit(9);
	 } */
	   
	 /* Accumulate values  */
	 if (Count<MAXTIP) {	   
	   isCal[Count] = rec[inOTF->myDesc->iloccal]!=0.0; /* Cal on? */
	   /* Get elevation in rad */
	   elev = ObitOTFArrayGeomElev (inOTF->geom, rec[inOTF->myDesc->iloct], 
					rec[inOTF->myDesc->ilocra], 
					rec[inOTF->myDesc->ilocdec]) * 1.745329e-2;
	   airmass[Count] = 1.0 / cos (1.5708-MAX(elev,0.1));  /* To airmass */
	   if (elev>minElR) {
	     for (j=0; j<ndetect; j++) 
	       data[j*MAXTIP+Count] = rec[inOTF->myDesc->ilocdata+j*incdatawt];
	   } else {
	     for (j=0; j<ndetect; j++) data[j*MAXTIP+Count] = fblank;
	   }
	   Count++;
	 } /* end  accumulate */
	 rec += inOTF->myDesc->lrec; /* Data record pointer */	 
       } /* end loop over buffer load */
  } /* end loop reading data */ 
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Better have some data */
  if (Count<=0) {
    Obit_log_error(err, OBIT_Error, "%s: NO data selected in %s", 
		   routine, inOTF->name);
    return;
 }

  /* Allocate TCal */
  TCal    = g_malloc0(ndetect*sizeof(float));

  /* Want plotting? */
  doPlot = FALSE;
  ObitInfoListGetTest (inOTF->info, "doPlot", &type, dim,  (gpointer)&doPlot);

  if (doPlot) {  /* Get info needed */
    /* TSky - default aTemp */
    TSky = ATEMP;
    ObitInfoListGetTest (inOTF->info, "tSky", &type, dim,  (gpointer)&TSky);
    
    /* tCal - no default */
    ObitInfoListGet (inOTF->info, "tCal", &type, dim,  (gpointer)TCal, err);
    if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
  }
 
  /* Loop over detectors */
  for (j=0; j<ndetect; j++) {

    /* Copy to work array - ObitOTFGetSolnAvgCal clobbers input */
    for (i=0;  i<Count; i++) work[i] = data[j*MAXTIP+i];
    /* get average cal value */
    cal = ObitOTFGetSolnAvgCal (Count, isCal, work);

    /* Correct cal data - convert to cal units */
    ical = 1.0 / cal;
    for (i=0; i<Count; i++) {
       if (data[j*MAXTIP+i]==fblank) continue;
       if (isCal[i]) data[j*MAXTIP+i] -= cal;
       data[j*MAXTIP+i] *= ical;   /* to cal units */
    }

    /* Fit scan */
    RMS[j] =  FitTip (Count, airmass, &data[j*MAXTIP], &ATemp[j], &TRx[j]);
  
    /* debug */
    fprintf (stdout, "detector %d ATemp %f TRx %f TSky %f TCal %f \n",j, 
	     ATemp[j], TRx[j], TSky, TCal[j]);

    /* Plot if requsted */
    if (doPlot && (RMS[j]!=fblank)) {
      /* Compute Tau */
      tau0 = TCal[j] * ATemp[j] / TSky;
      PlotTip (Count, scan, j+1, tau0, airmass, &data[j*MAXTIP], ATemp[j], TRx[j], 
	       err);
      if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
    } /* end plotting */
    
  } /* end loop over detectors */

  /* Cleanup */
  if (TCal)    g_free(TCal);
  if (data)    g_free(data);
}  /* end FitTipScan */

/**
 * Fits tipping scan to determine receiver temperature and sky brightness
 * Assumes sky brightness increases linearly with air mass (breaks down 
 * for large opacity).
 *  T = Trx + Tatm * airmass, where airmass = 1.0 / cos (za)       
 * \param ndata    number of data points 
 * \param airmass  Array of airmass
 * \param Temp     Array of data points (cal units)
 * \param Tatm    [out] The sky brightness in cal units per airmass
 * \param Trx     [out] Receiver temperature in units of the cal
 * \return RMS residual (cal units), blanked if undetermined.
 */
  static float 
    FitTip (olong ndata, ofloat *airmass, ofloat *Temp, ofloat *Tatm, ofloat *Trx)
{
  odouble x, y, r, Sxx, Sxy;
  odouble sumxT=0.0, sumyT=0.0, sumxxT=0.0, sumyyT=0.0, sumxyT=0.0;
  ofloat rms, fblank=ObitMagicF();
  olong i, countT=0;

  /* Do sums */
  for (i=0; i<ndata; i++) {
    x = airmass[i];
    if (Temp[i]!=fblank) {
      y = Temp[i];
      countT++;
      sumxT += x;
      sumyT += y;
      sumxxT += x*x;
      sumyyT += y*y;
      sumxyT += x*y;
    }
  } /* end loop over data */
  
  /* Determine coefficients */
  if (countT>0) {
    Sxx = sumxxT - sumxT*sumxT/countT;
    Sxy = sumxyT - sumxT*sumyT/countT;
    *Tatm = Sxy/Sxx;
    *Trx  = (sumyT-(*Tatm)*sumxT)/countT;
  } else {
    *Tatm = fblank;
    *Trx  = fblank;
  }

  /* Get RMS residual */
  countT = 0;
  sumxxT = 0.0;
  for (i=0; i<ndata; i++) {
    if (Temp[i]!=fblank) {
      r = Temp[i] - (*Trx) - (*Tatm)*airmass[i];
      countT++;
      sumxxT += r*r;
    }
  } /* end loop over data */
  if (countT > 0) rms = (ofloat)sqrt(sumxxT/countT);
  else rms = fblank;

  return rms;
} /* end FitTip */

/**
 * Plot data plus model
 * \param Count  Number of samples to be plotted
 * \param Time   Time tags of the data 
 * \param scan   scan number
 * \param detect detector number (1-rel)
 * \param Temp   Measured values
 * \param base   value off source
 * \param peak   peak value
 * \param center Time of center of Gaussian
 * \param sigma  Gaussian sigma (width) in time units.
 * \param err    Error stack, returns if not empty.
 */
static void PlotCalGauss (olong Count, ofloat *Time,  
			  olong scan, olong detect, ofloat *Temp,
			  ofloat base, ofloat peak, ofloat center, 
			  ofloat sigma, ObitErr *err)
{
  ofloat *pdata=NULL, *px=NULL, ivar, ymax, ymin, fblank = ObitMagicF();
  olong i, nplot;
  ObitPlot *plot = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar strTemp[121];
  gchar *routine = "PlotCalGauss";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  if (Count<=0) return;

  /* initialize plotting */
  plot = newObitPlot (routine);

  /* Plot labels */
  g_snprintf (strTemp, 120, "Data and model scan %d detector %d",scan, detect);
  dim[0] = strlen (strTemp);
  ObitInfoListAlwaysPut (plot->info, "TITLE", OBIT_string, dim, strTemp);
  strncpy (strTemp, "Time (sec)", 120);
  dim[0] = strlen (strTemp);
  ObitInfoListAlwaysPut (plot->info, "XLABEL", OBIT_string, dim, strTemp);
  strncpy (strTemp, "Offset", 120);
  dim[0] = strlen (strTemp);
  ObitInfoListAlwaysPut (plot->info, "YLABEL", OBIT_string, dim, strTemp);

  /* initialize plotting */
  ObitPlotInitPlot (plot, NULL, 15, 1,1, err);
  if (err->error) {
    plot = ObitPlotUnref(plot);
    Obit_traceback_msg (err, routine, plot->name);
  }

  /* Find extrema */
  ymax = -1.0e20;
  ymin =  1.0e20;
  for (i=0; i<Count; i++) {
    if (Temp[i]!=fblank) {
      ymax = MAX (ymax, Temp[i]);
      ymin = MIN (ymin, Temp[i]);
    }
  }
  /* Set extrema */
  dim[0] = 1;
  ObitInfoListAlwaysPut(plot->info, "YMAX", OBIT_float, dim, (gpointer*)&ymax);
  ObitInfoListAlwaysPut(plot->info, "YMIN", OBIT_float, dim, (gpointer*)&ymin);

  /* Allocate model storage */
  pdata = g_malloc0(Count*sizeof(float));
  px    = g_malloc0(Count*sizeof(float));
  
  /* model to plot */
  nplot = 0;
  if (sigma!=0.0) ivar = 1.0 / (sigma*sigma); /* inverse variance of Gaussian */
  else ivar = 1.0;
  for (i=0; i<Count; i++) {
    pdata[nplot] = base + peak * exp (-(Time[i]-center)*(Time[i]-center)*ivar);
    px[nplot] = Time[i];
    nplot++;
  }
  /* plot it */
  ObitPlotXYPlot (plot, -1, nplot, px, pdata, err);

  /* get  valid data to plot */
  nplot = 0;
  for (i=0; i<Count; i++) {
    if (Temp[i]!=fblank) {
      pdata[nplot++] = Temp[i];
      px[nplot] = Time[i];
    }
  }

  /* plot it */
  ObitPlotXYOver (plot, 2, nplot, px, pdata, err);
  ObitPlotFinishPlot (plot, err);  /* Finalize plot */

  /* cleanup */
  plot = ObitPlotUnref(plot);
  if (pdata) g_free(pdata);  /* Deallocate plot memory */
  if (px)    g_free(px); 

  if (err->error) Obit_traceback_msg (err, routine, plot->name);
} /* end PlotCalGauss */

/**
 * Plot tip data plus model
 *  T = Trx + Tatm * airmass, where airmass = 1.0 / cos (za)       
 * \param Count    Number of samples to be plotted
 * \param scan     Scan number
 * \param detect   Detector number (1-rel)
 * \param tau0     Estimated atmospheric opacity
 * \param airmass  Array of airmass
 * \param Temp     Array of data points (cal units)
 * \param Tatm     The sky brightness in cal units per airmass
 * \param Trx      Receiver temperature in units of the cal
 * \param err      Error stack, returns if not empty.
 */
static void PlotTip (olong Count, olong scan, olong detect, ofloat tau0,
		     ofloat *airmass, ofloat *Temp, ofloat Tatm, 
		     ofloat Trx, ObitErr *err)
{
  ofloat *pdata=NULL, *px=NULL, ymax, ymin, fblank = ObitMagicF();
  olong i, nplot;
  ObitPlot *plot = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar strTemp[121];
  gchar *routine = "PlotTip";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  if (Count<=0) return;

  /* initialize plotting */
  plot = newObitPlot (routine);

  /* Plot labels */
  g_snprintf (strTemp, 120, "Sky Tip scan %d detector %d tau0=%f",
	      scan, detect, tau0);
  dim[0] = strlen (strTemp);
  ObitInfoListAlwaysPut (plot->info, "TITLE", OBIT_string, dim, strTemp);
  strncpy (strTemp, "Air Mass", 120);
  dim[0] = strlen (strTemp);
  ObitInfoListAlwaysPut (plot->info, "XLABEL", OBIT_string, dim, strTemp);
  strncpy (strTemp, "Offset in cal units", 120);
  dim[0] = strlen (strTemp);
  ObitInfoListAlwaysPut (plot->info, "YLABEL", OBIT_string, dim, strTemp);

  /* initialize plotting */
  ObitPlotInitPlot (plot, NULL, 15, 1,1, err);
  if (err->error) {
    plot = ObitPlotUnref(plot);
    Obit_traceback_msg (err, routine, plot->name);
  }

  /* Find extrema */
  ymax = -1.0e20;
  ymin =  1.0e20;
  for (i=0; i<Count; i++) {
    if (Temp[i]!=fblank) {
      ymax = MAX (ymax, Temp[i]);
      ymin = MIN (ymin, Temp[i]);
    }
  }
  /* Set extrema */
  dim[0] = 1;
  ObitInfoListAlwaysPut(plot->info, "YMAX", OBIT_float, dim, (gpointer*)&ymax);
  ObitInfoListAlwaysPut(plot->info, "YMIN", OBIT_float, dim, (gpointer*)&ymin);

  /* Allocate model storage */
  pdata = g_malloc0(Count*sizeof(float));
  px    = g_malloc0(Count*sizeof(float));
  
  /* model to plot */
  nplot = 0;
  for (i=0; i<Count; i++) {
    pdata[nplot] = Trx + Tatm * airmass[i];
    px[nplot] = airmass[i];
    nplot++;
  }
  /* plot it */
  ObitPlotXYPlot (plot, -1, nplot, px, pdata, err);

  /* get  valid data to plot */
  nplot = 0;
  for (i=0; i<Count; i++) {
    if (Temp[i]!=fblank) {
      pdata[nplot++] = Temp[i];
      px[nplot] = airmass[i];
    }
  }

  /* plot it */
  ObitPlotXYOver (plot, 2, nplot, px, pdata, err);
  ObitPlotFinishPlot (plot, err);  /* Finalize plot */

  /* cleanup */
  plot = ObitPlotUnref(plot);
  if (pdata) g_free(pdata);  /* Deallocate plot memory */
  if (px)    g_free(px); 

  if (err->error) Obit_traceback_msg (err, routine, plot->name);
} /* end PlotTip */


