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

#include "ObitOTFGetAtmCor.h"
#include "ObitOTFGetSoln.h"
#include "ObitPennArrayAtmFit.h"
#include "ObitGBTCCBUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFGetAtmCor.c
 * ObitOTFGetAtmCor utility module calibration function definitions.
 */

/*---------------Private function prototypes----------------*/
/** return opacity */
static ofloat getTau0 (ofloat time, gint32 taudim[], ofloat *taudata);

/*----------------------Public functions---------------------------*/
/**
 * Determine gain offset calibration for an OTF dataset.
 * The offset is determined from an atmospheric model.
 * The gain calibration is determined from the average noise cal values
 * and the atmospheric opacity.
 * Calibration parameters are on the inOTF info member.
 * \li "solInt"   OBIT_float (1,1,1) Solution interval in days [def 10 sec].
 * \li "Tau0"     OBIT_float (?,?,1) Zenith opacity in nepers [def 0].
 *                Can be passed as either a constant scalar or
 *                an array of (time, tau0) pairs.
 * \li "aTemp"    OBIT_float (*,1,1) Effective atmospheric temperature in data units.
 *                i.e. the additional offset per airmass due to the atmosphere.
 *                per detector in units of the cal.  
 *                If only one given, it is used for all detectors
 * \li "minEl"    OBIT_float (1,1,1) Minimum elevation (deg)
 * \li "tRx"      OBIT_float (*,1,1) Receiver temperature per detector in units of the cal
 *                If only one given, it is used for all detectors
 * \li "calJy"    OBIT_float (*,1,1) Noise cal value in Jy, per detector [def 1.0] .
 *                If only one given, it is used for all detectors
 * \li "Azoff"    OBIT_float (*,1,1) Constant offset in deg to add to Az [def 0] .
 * \li "Eloff"    OBIT_float (*,1,1) Constant offset in deg to add to El  [def 0.0] .
 * \param inOTF    Input OTF data. 
 * \param outOTF   OTF with which the output  OTFSoln is to be associated
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created OTFSoln object which is associated with inOTF.
 */
ObitTableOTFSoln* ObitOTFGetAtmCor (ObitOTF *inOTF, ObitOTF *outOTF, ObitErr *err)
{
#define MAXSAMPLE 1000   /* Maximum number of samples in an integration */
  ObitTableOTFSoln *outSoln=NULL;
  ObitTableOTFSolnRow *row=NULL;
  gint32 taudim[MAXINFOELEMDIM] = {1,1,1,1,1}, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat *rec, solInt, *tau0, minEl, t0, tend, el=0.0, el0=0.0, elend=0.0, tMult, airmass;
  ofloat Azoff, Eloff, lastScan, lastTarget;
  ofloat *atemp=NULL, *trx=NULL, *calJy=NULL;
  ofloat *calArray=NULL, calFlag=0.0, fblank = ObitMagicF();
  olong iRow, ver, i, j, lrec, ncoef;
  gboolean flag, someData=FALSE, *isCal=NULL, doCalSelect;
  olong  npoly, ndetect, calCount, incdatawt;
  ObitIOCode retCode;
  ObitIOAccess access;
  ObitOTFArrayGeom* geom;
  gchar *tname;
  gchar *routine = "ObitOTFGetAtmCor";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outSoln;
  g_assert (ObitOTFIsA(inOTF));
  g_assert (ObitOTFArrayGeomIsA(inOTF->geom));

  /* open OTF data to fully instantiate if not already open */
  ObitInfoListGetTest(inOTF->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;
  if ((inOTF->myStatus==OBIT_Inactive) || (inOTF->myStatus==OBIT_Defined)) {
    retCode = ObitOTFOpen (inOTF, access, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
  }

  /* Get array geometry */
  geom = inOTF->geom;

  /* Create output */
  tname = g_strconcat ("Calibration for: ",inOTF->name, NULL);
  ver = 0;
  ncoef = 1; 
  npoly = ncoef;
  outSoln = newObitTableOTFSolnValue(tname, (ObitData*)outOTF, &ver, 
				     OBIT_IO_WriteOnly,  
				     inOTF->geom->numberDetect, npoly, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Create work arrays */
  lrec = inOTF->myDesc->lrec;
  t0 = -1.0e20;
  tend = t0;
  elend = 45.0;
  lastScan = -1000.0;
  lastTarget = 0.0;
  ndetect = inOTF->geom->numberDetect;  /* number of detectors */
  trx   = g_malloc0(ndetect*sizeof(ofloat));
  atemp = g_malloc0(ndetect*sizeof(ofloat));
  calJy = g_malloc0(ndetect*sizeof(ofloat));
  calArray = g_malloc0(ndetect*MAXSAMPLE*sizeof(ofloat));
  isCal    = g_malloc0(MAXSAMPLE*sizeof(gboolean));

  /* Get parameters for calibration */
  /* Solution interval default 10 sec */
  solInt = 10.0 / 86400.0;
  ObitInfoListGetTest(inOTF->info, "solInt", &type, dim, (gpointer*)&solInt);
  tau0 = NULL;
  ObitInfoListGetP(inOTF->info, "Tau0",   &type, taudim, (gpointer*)&tau0);
  minEl = 1.0;
  ObitInfoListGetTest(inOTF->info, "minEl",  &type, dim, (gpointer*)&minEl);
  for (j=0; j<ndetect; j++) trx[j] = 1.0;
  ObitInfoListGetTest(inOTF->info, "tRx",  &type, dim, (gpointer*)trx);
  /* If only one given, use it for all */
  if (dim[0]==1) for (j=1; j<ndetect; j++) trx[j] = trx[0];
  for (j=0; j<ndetect; j++) atemp[j] = 300.0;
  ObitInfoListGetTest(inOTF->info, "aTemp",  &type, dim, (gpointer*)atemp);
  /* If only one given, use it for all */
  if (dim[0]==1) for (j=1; j<ndetect; j++) atemp[j] = atemp[0];
  for (j=0; j<ndetect; j++) calJy[j] = 1.0;
  ObitInfoListGetTest(inOTF->info, "calJy",  &type, dim, (gpointer*)calJy);
  if (dim[0]==1) for (j=1; j<ndetect; j++) calJy[j] = calJy[0];
  Azoff = 0.0;
  ObitInfoListGetTest(inOTF->info, "Azoff",  &type, dim, (gpointer*)&Azoff);
  Eloff = 0.0;
  ObitInfoListGetTest(inOTF->info, "Eloff", &type, dim, (gpointer*)&Eloff);
  incdatawt = inOTF->myDesc->incdatawt; /* increment in data-wt axis */
 
  /* Open table */
  if ((ObitTableOTFSolnOpen (outSoln, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output OTFSoln table");
    return outSoln;
  }

  /* Create Row */
  row = newObitTableOTFSolnRow (outSoln);

  /* Attach row to output buffer */
  ObitTableOTFSolnSetRow (outSoln, row, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Initialize */
  row->dAz =  Azoff;  /* any pointing corrections */
  row->dEl =  Eloff;  /* any pointing corrections */
  row->Target = 0;
  for (j=0; j<npoly; j++)   row->poly[j] = 0.0;
  for (j=0; j<ndetect; j++) row->mult[j] = 1.0;
  for (j=0; j<ndetect; j++) row->wt[j]   = 1.0;
  for (j=0; j<ndetect; j++) row->add[j]  = 0.0;
  flag = FALSE;

  /* init cal average arrays */
  calCount = 0;      /* Cal values in array */

  /* loop calibrating data */
  retCode = OBIT_IO_OK;
  while (retCode == OBIT_IO_OK) {

    /* read buffer */
    if (doCalSelect) retCode = ObitOTFReadSelect (inOTF, NULL, err);
    else retCode = ObitOTFRead (inOTF, NULL, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
    if (retCode==OBIT_IO_EOF) break; /* done? */

    /* Record pointer */
    rec = inOTF->buffer;
  
    /* Loop over buffer */
       for (i=0; i<inOTF->myDesc->numRecBuff; i++) {

	 /* elevation */
	 el =  ObitOTFArrayGeomElev (geom, rec[inOTF->myDesc->iloct], 
				    rec[inOTF->myDesc->ilocra], 
				    rec[inOTF->myDesc->ilocdec]);
	 /* Need to set first time/elevation in integration? */
	 if (t0<-1.0e10) {
	   t0 = rec[inOTF->myDesc->iloct];
	   tend = t0;
	   el0 =  el;
	   lastScan = rec[inOTF->myDesc->ilocscan];  /* Which scan number */
	   lastTarget = rec[inOTF->myDesc->iloctar]; /* Which target number */
	 }

	 /* Time for a new correction? If so, compute correction and write.*/
	 if ((rec[inOTF->myDesc->iloct] > (t0+solInt)) || 
	     /* or end of scan */
	     ( rec[inOTF->myDesc->ilocscan] != lastScan) ||
	     /* Or buffers full */
	     (calCount>=MAXSAMPLE) || 
	     /* or end of data */
	     /* wrong??? (i+1 >= inOTF->myDesc->numRecBuff)) {*/
	     (i+1+inOTF->myDesc->firstRec >= inOTF->myDesc->nrecord)) {

	   /* Set descriptive info on Row */
	   row->Time  = t0 - 0.001/86400.0; /* minus msec */
	   row->TimeI = solInt;
	   row->Target = (oint)(lastTarget+0.5);

	   /* Average cal values */
	   for (j=0; j<ndetect; j++) {
	     row->cal[j] = ObitOTFGetSolnAvgCal (calCount, isCal, 
						 &calArray[j*MAXSAMPLE]);
	   }

	   /* If this is the CCB some massaging of the cals is required */
	   if (inOTF->myDesc->OTFType==OBIT_GBTOTF_CCB) 
	     ObitGBTCCBAvgCal(inOTF->myDesc, calFlag, ndetect, row->cal);

	   /* record at beginning of solution interval */
	   /* Number of atmospheres seen through */
	   if (el0 > 0.1) 
	     airmass = 1.0 / cos (DG2RAD * (90.0 - el0));
	   else
	     airmass = 10.0;

	   tMult = exp (getTau0(row->Time,taudim,tau0) * airmass); /* multiplicative, additive terms */

           /* set calibration part of Soln */
	   for (j=0; j<ndetect; j++) {
	     if (row->cal[j] != fblank) {
	       row->mult[j] = tMult * calJy[j] / row->cal[j];
	       row->add[j]  = -(atemp[j] * airmass + trx[j]) * row->cal[j];
	     } else {
	       row->mult[j] = fblank;
	       row->add[j]  = fblank;
	     }
	   }

	   /* flag solution if too low elevation */
	   if (flag) {
	     for (j=0; j<ndetect; j++) row->mult[j] = fblank;
	     for (j=0; j<ndetect; j++) row->add[j]  = fblank;
	   }

	   /* Write Soln table */
	   iRow = -1;
	   if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
		!= OBIT_IO_OK) || (err->error>0)) { 
	     Obit_log_error(err, OBIT_Error, "ERROR writing OTFSoln Table file");
	     return outSoln;
	   }

	   /* record at end of solution interval */
	   row->Time  = tend + 0.001/86400.0; /* time plus msec */
	   /* Number of atmospheres seen through */
	   if (elend > 0.1) 
	     airmass = 1.0 / cos (DG2RAD * (90.0 - elend));
	   else
	     airmass = 10.0;

	   tMult = exp (getTau0(row->Time,taudim,tau0) * airmass); /* multiplicative, additive terms */

           /* set calibration part of Soln */
	   for (j=0; j<ndetect; j++) {
	     if (row->cal[j] != fblank) {
	       row->mult[j] = tMult * calJy[j] / row->cal[j];
	       row->add[j]  = -(atemp[j] * airmass + trx[j]) * row->cal[j];
	     } else {
	       row->mult[j] = fblank;
	       row->add[j]  = fblank;
	     }
	   }

	   /* flag solution if too low elevation */
	   if (flag) {
	     for (j=0; j<ndetect; j++) row->mult[j] = fblank;
	     for (j=0; j<ndetect; j++) row->add[j]  = fblank;
	   }
	   flag = FALSE;  /* reset flag */

	   /* Write Soln table */
	   iRow = -1;
	   if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
		!= OBIT_IO_OK) || (err->error>0)) { 
	     Obit_log_error(err, OBIT_Error, "ERROR writing OTFSoln Table file");
	     return outSoln;
	   }

	   /* start new solution interval */
	   t0 = rec[inOTF->myDesc->iloct];
	   el0 = el;

	   /* reset cal sums */
	   calCount = 0; /* Cal values in sumArray */

	   /* Done? */
	   if (i+1+inOTF->myDesc->firstRec >= inOTF->myDesc->nrecord) break;
	 } /* end do solution */
	 
	 /* Accumulate values for determining the cal value */
	 if (calCount<MAXSAMPLE) {
	   isCal[calCount] = rec[inOTF->myDesc->iloccal]!=0.0; /* Cal on? */
	   /* What is the value of the cal flag? */
	   if ((isCal[calCount]) && (calFlag==0.0)) 
	     calFlag = rec[inOTF->myDesc->iloccal];
	   for (j=0; j<ndetect; j++) 
	     calArray[j*MAXSAMPLE+calCount] = 
	       rec[inOTF->myDesc->ilocdata+j*incdatawt];
	   calCount++;
	   someData = TRUE;
	 } /* end save values for cal average */

	 /* Is some data below the limit? */
	 flag = flag || (el < minEl);

	 lastScan = rec[inOTF->myDesc->ilocscan];  /* Which scan number */
	 lastTarget = rec[inOTF->myDesc->iloctar]; /* Which target number */
	 elend = el;                       /* last elevation in solution */
	 tend = rec[inOTF->myDesc->iloct]; /* last time in solution */
	 rec += inOTF->myDesc->lrec;       /* Data record pointer */
	 
       } /* end loop over buffer load */
  } /* end loop reading data */
  
  /* Anything left in accumulators?  */
  if (calCount>0) {
    /* Set descriptive info on Row */
    row->Time  = t0 - 0.001/86400.0; /* minus msec */
    row->TimeI = solInt;
    row->Target = (oint)(lastTarget+0.5);
    
    /* Average cal values */
    for (j=0; j<ndetect; j++) {
      row->cal[j] = ObitOTFGetSolnAvgCal (calCount, isCal, &calArray[j*MAXSAMPLE]);
    }
    
    /* If this is the CCB some massaging of the cals is required */
    if (inOTF->myDesc->OTFType==OBIT_GBTOTF_CCB) 
      ObitGBTCCBAvgCal(inOTF->myDesc, calFlag, ndetect, row->cal);
    
    /* record at beginning of solution interval */
    /* Number of atmospheres seen through */
    if (el0 > 0.1) 
      airmass = 1.0 / cos (DG2RAD * (90.0 - el0));
    else
      airmass = 10.0;
    
    tMult = exp (getTau0(row->Time,taudim,tau0) * airmass); /* multiplicative, additive terms */
    
    /* set calibration part of Soln */
    for (j=0; j<ndetect; j++) {
      if (row->cal[j] != fblank) {
	row->mult[j] = tMult * calJy[j] / row->cal[j];
	row->add[j]  = -(atemp[j] * airmass + trx[j]) * row->cal[j];
      } else {
	row->mult[j] = fblank;
	row->add[j]  = fblank;
      }
    }
    
    /* flag solution if too low elevation */
    if (flag) {
      for (j=0; j<ndetect; j++) row->mult[j] = fblank;
      for (j=0; j<ndetect; j++) row->add[j]  = fblank;
    }
    
    /* Write Soln table */
    iRow = -1;
    if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR writing OTFSoln Table file");
      return outSoln;
    }
    
    /* record at end of solution interval */
    row->Time  = tend + 0.001/86400.0; /* time plus msec */
    /* Number of atmospheres seen through */
    if (elend > 0.1) 
      airmass = 1.0 / cos (DG2RAD * (90.0 - elend));
    else
      airmass = 10.0;
    
    tMult = exp (getTau0(row->Time,taudim,tau0) * airmass);  /* multiplicative, additive terms */
    
    /* set calibration part of Soln */
    for (j=0; j<ndetect; j++) {
      if (row->cal[j] != fblank) {
	row->mult[j] = tMult * calJy[j] / row->cal[j];
	row->add[j]  = -(atemp[j] * airmass + trx[j]) * row->cal[j];
      } else {
	row->mult[j] = fblank;
	row->add[j]  = fblank;
      }
    }
    
    /* flag solution if too low elevation */
    if (flag) {
      for (j=0; j<ndetect; j++) row->mult[j] = fblank;
      for (j=0; j<ndetect; j++) row->add[j]  = fblank;
    }
    flag = FALSE;  /* reset flag */
    
    /* Write Soln table */
    iRow = -1;
    if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR writing OTFSoln Table file");
      return outSoln;
    }
    
    /* start new solution interval */
    t0 = rec[inOTF->myDesc->iloct];
    el0 = el;
    
  } /* end final write */

  /* Close cal table */
  if ((ObitTableOTFSolnClose (outSoln, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output OTFSoln Table file");
    return outSoln;
  }
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Give warning if no data selected */
  if (!someData) Obit_log_error(err, OBIT_InfoWarn, 
				"%s: Warning: NO data selected", routine);

  /* Cleanup */
  if (trx) g_free(trx);
  if (atemp) g_free(atemp);
  if (calJy) g_free(calJy);
  if (calArray) g_free(calArray);
  if (isCal) g_free(isCal);

  return outSoln;
} /* end ObitOTFGetAtmCor */

/**
 * Create OTFSoln table putting atmospheric emission corrections into the detector offsets.
 * \li "solInt"    OBIT_float (1,1,1) Solution interval in days [def 1 sec].
 * \li "AtmEm"     OBIT_float (1,1,1) Equivalent zenith atmospheric brightness in Jy [0]
 * \li "Tau0"      OBIT_float (?,?,1) Zenith opacity in nepers [def 0].
 *                 Can be passed as either a constant scalar or
 *                 an array of (time, tau0) pairs.
 * \param inOTF    Input OTF data. 
 * \param outOTF   OTF with which the output  OTFSoln is to be associated
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created OTFSoln object which is associated with outOTF.
 */
ObitTableOTFSoln* ObitOTFGetAtmEm (ObitOTF *inOTF, ObitOTF *outOTF, ObitErr *err)
{
  ObitTableOTFSoln *outSoln=NULL;
  ObitTableOTFSolnRow *row=NULL;
  ObitOTFDesc *desc=NULL;
  ObitIOAccess access;
  gint32 taudim[MAXINFOELEMDIM] = {1,1,1,1,1}, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat *rec, solInt, t0, sumTime, lastTime=-1.0, lastTarget=-1.0, lastScan=-1.0;
  ofloat *tau0, tMult, tAdd, AtmEm, el, airmass, sumRA=0.0, sumDec=0.0, lastRA=0.0, lastDec=0.0;
  ofloat ra0=-1.0e20, dec0=-1.0e20, tAdd0;
  olong iRow, i, j, lrec, ver, ncoef;
  olong  nDet, nTime, npoly;
  gboolean doCalSelect, someData=FALSE;
  ObitIOCode retCode;
  gchar *tname;
  gchar *routine = "ObitOTFGetAtmEm";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outSoln;
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
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
  }
  lrec = inOTF->myDesc->lrec;
  nDet = inOTF->geom->numberDetect;
  t0 = -1.0e20;

  /* Create output */
  tname = g_strconcat ("Calibration for: ",inOTF->name, NULL);
  ver = 0;
  ncoef = 1;  /* do first order atmospheric model (0=1, 1=3, 2nd = 6, 3rd=10) */
  npoly = ncoef;
  outSoln = newObitTableOTFSolnValue(tname, (ObitData*)outOTF, &ver, OBIT_IO_WriteOnly,  
				   inOTF->geom->numberDetect, npoly, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Get parameters for calibration */
  /* "Solution interval" default 1 sec */
  solInt = 1.0 / 86400.0;
  ObitInfoListGetTest(inOTF->info, "solInt", &type, dim, &solInt);

  /* Opacity */
  tau0 = NULL;
  ObitInfoListGetP(inOTF->info, "Tau0",   &type, taudim, (gpointer*)&tau0);

  /* Zenith Emission */
  AtmEm = 0.0;
  ObitInfoListGetTest(inOTF->info, "AtmEm",   &type, dim, &AtmEm);

  /* Open table */
  if ((ObitTableOTFSolnOpen (outSoln, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening input OTFSoln table", routine);
    return outSoln;
  }

  /* Create Row */
  row = newObitTableOTFSolnRow (outSoln);

  /* Attach row to output buffer */
  ObitTableOTFSolnSetRow (outSoln, row, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Initialize */
  row->dAz = 0.0;  /* no pointing corrections */
  row->dEl = 0.0;  /* no pointing corrections */
  for (i=0; i<nDet; i++) {
    row->mult[i] = 1.0;
    row->wt[i]   = 1.0;
    row->add[i]  = 0.0;
    row->cal[i]  = 0.0;
  }

  /* loop calibrating data */
  retCode = OBIT_IO_OK;
  sumTime = 0.0;
  sumRA   = 0.0;
  sumDec  = 0.0;
  nTime = 0;
  while (retCode == OBIT_IO_OK) {
    
    /* read buffer */
    if (doCalSelect) retCode = ObitOTFReadSelect (inOTF, NULL, err);
    else retCode = ObitOTFRead (inOTF, NULL, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);
    if (retCode==OBIT_IO_EOF) break; /* done? */
    
    /* Record pointer */
    rec = inOTF->buffer;
    
    /* First time */
    if (t0<-1.0e10) {
      t0         = rec[inOTF->myDesc->iloct];
      lastTarget = rec[inOTF->myDesc->iloctar];
      lastScan   = rec[inOTF->myDesc->ilocscan];
      lastTime   = rec[inOTF->myDesc->iloct];
      lastRA     = rec[inOTF->myDesc->ilocra];
      lastDec    = rec[inOTF->myDesc->ilocdec];
      if (ra0<1000.0) {
	ra0  = lastRA;
	dec0 = lastDec;
      }
      /* Write entry at very beginning */
      row->Time  = rec[inOTF->myDesc->iloct]; 
      row->TimeI = 0.0;
      row->Target = (oint)(rec[inOTF->myDesc->iloctar]+0.5);
      iRow = -1;

      /*  differential Atmospheric correction */
      el = ObitOTFArrayGeomElev (inOTF->geom, row->Time, ra0, dec0);
      if (el > 0.1) 
	airmass = 1.0 / cos (DG2RAD * (90.0 - el));
      else
	airmass = 10.0;
      tMult = exp (-getTau0(row->Time,taudim,tau0) * airmass); /* multiplicative term */
      tAdd0  = AtmEm * (1.0 - tMult);      /* additive term */
      el = ObitOTFArrayGeomElev (inOTF->geom, row->Time, 
				 rec[inOTF->myDesc->ilocra], 
				 rec[inOTF->myDesc->ilocdec]);
      if (el > 0.1) 
	airmass = 1.0 / cos (DG2RAD * (90.0 - el));
      else
	airmass = 10.0;
      tMult = exp (-getTau0(row->Time,taudim,tau0) * airmass); /* multiplicative term */
      tAdd  = AtmEm * (1.0 - tMult);      /* additive term */
      for (j=0; j<nDet; j++) row->add[j] = -tAdd + tAdd0;
	
      if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "%s: ERROR writing OTFSoln Table file", routine);
	return outSoln;
      }
    }
    
    /* Loop over buffer */
    for (i=0; i<inOTF->myDesc->numRecBuff; i++) {

      /* Accumulation or scan finished? If so, write "calibration".*/
      if ((rec[inOTF->myDesc->iloct] > (t0+solInt)) || 
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

	    /*  differential Atmospheric correction */
	    el = ObitOTFArrayGeomElev (inOTF->geom, row->Time, ra0, dec0);
	    if (el > 0.1) 
	      airmass = 1.0 / cos (DG2RAD * (90.0 - el));
	    else
	      airmass = 10.0;
	    tMult = exp (-getTau0(row->Time,taudim,tau0) * airmass); /* multiplicative term */
	    tAdd0  = AtmEm * (1.0 - tMult);      /* additive term */

	    el = ObitOTFArrayGeomElev (inOTF->geom, row->Time, lastRA, lastDec);
	    if (el > 0.1) 
	      airmass = 1.0 / cos (DG2RAD * (90.0 - el));
	    else
	      airmass = 10.0;
	    tMult = exp (-getTau0(row->Time,taudim,tau0) * airmass); /* multiplicative term */
	    tAdd  = AtmEm * (1.0 - tMult);      /* additive term */
	    for (j=0; j<nDet; j++) row->add[j] = -tAdd + tAdd0;
	
	    if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
		 != OBIT_IO_OK) || (err->error>0)) { 
	      Obit_log_error(err, OBIT_Error, "%s: ERROR writing OTFSoln Table file", routine);
	      return outSoln;
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
      
	    /*  differential Atmospheric correction */
	    el = ObitOTFArrayGeomElev (inOTF->geom, row->Time, ra0, dec0);
	    if (el > 0.1) 
	      airmass = 1.0 / cos (DG2RAD * (90.0 - el));
	    else
	      airmass = 10.0;
	    tMult = exp (-getTau0(row->Time,taudim,tau0) * airmass); /* multiplicative term */
	    tAdd0  = AtmEm * (1.0 - tMult);      /* additive term */

	    el = ObitOTFArrayGeomElev (inOTF->geom, row->Time, sumRA/nTime, sumDec/nTime);
	    if (el > 0.1) 
	      airmass = 1.0 / cos (DG2RAD * (90.0 - el));
	    else
	      airmass = 10.0;
	    tMult = exp (-getTau0(row->Time,taudim,tau0) * airmass); /* multiplicative term */
	    tAdd  = AtmEm * (1.0 - tMult);      /* additive term */
	    for (j=0; j<nDet; j++) row->add[j] = -tAdd + tAdd0;
	  }

	  /* Write Soln table */
	  iRow = -1;
	  if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	       != OBIT_IO_OK) || (err->error>0)) { 
	    Obit_log_error(err, OBIT_Error, "%s: ERROR writing OTFSoln Table file", routine);
	    return outSoln;
	  }
	  /* initialize accumulators */
	  t0         = rec[inOTF->myDesc->iloct];
	  lastTarget = rec[inOTF->myDesc->iloctar];
	  lastScan   = rec[inOTF->myDesc->ilocscan];
	  lastRA     = rec[inOTF->myDesc->ilocra];
	  lastDec    = rec[inOTF->myDesc->ilocdec];
	  sumTime    = 0.0;
	  sumRA      = 0.0;
	  sumDec     = 0.0;
	  nTime      = 0;
	} /* end of do solution if there is data */
	
      } /* end do solution */
      
      /* accumulate */
      lastTime = rec[inOTF->myDesc->iloct];
      lastRA   = rec[inOTF->myDesc->ilocra];
      lastDec  = rec[inOTF->myDesc->ilocdec];
      nTime++; /* how many data points */
      sumTime += rec[inOTF->myDesc->iloct];
      sumRA    += rec[inOTF->myDesc->ilocra];
      sumDec   += rec[inOTF->myDesc->ilocdec];
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
    /* Write Soln table */
    iRow = -1;

    /*  differential Atmospheric correction */
    el = ObitOTFArrayGeomElev (inOTF->geom, row->Time, ra0, dec0);
    if (el > 0.1) 
      airmass = 1.0 / cos (DG2RAD * (90.0 - el));
    else
      airmass = 10.0;
    tMult = exp (-getTau0(row->Time,taudim,tau0) * airmass); /* multiplicative term */
    tAdd0  = AtmEm * (1.0 - tMult);      /* additive term */

    el = ObitOTFArrayGeomElev (inOTF->geom, row->Time, lastRA, lastDec);
    if (el > 0.1) 
      airmass = 1.0 / cos (DG2RAD * (90.0 - el));
    else
      airmass = 10.0;
    tMult = exp (-getTau0(row->Time,taudim,tau0) * airmass); /* multiplicative term */
    tAdd  = AtmEm * (1.0 - tMult);      /* additive term */
    for (j=0; j<nDet; j++) row->add[j] = -tAdd + tAdd0;
	
    if ((ObitTableOTFSolnWriteRow (outSoln, iRow, row, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "%s: ERROR writing OTFSoln Table file", routine);
      return outSoln;
    }
  } /* End final cal */

  /* Close soln table */
  if ((ObitTableOTFSolnClose (outSoln, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing input OTFSoln Table file", routine);
    return outSoln;
  }
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outSoln);

  /* Give warning if no data selected */
  if (!someData) Obit_log_error(err, OBIT_InfoWarn, 
				"%s: Warning: NO data selected", routine);

  return outSoln;
} /* end ObitOTFGetAtmEm */

/*----------------------Private functions---------------------------*/

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
