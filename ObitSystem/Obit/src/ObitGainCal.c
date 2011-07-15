/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010                                               */
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
#include "ObitGainCal.h"
#include "ObitUVDesc.h"
#include "ObitTableNX.h"
#include "ObitTableCL.h"
#include "ObitTableSN.h"
#include "ObitTableSUUtil.h"
#include "ObitTableANUtil.h"
#include "ObitWeather.h"
#include "ObitOpacity.h"
#include "ObitTsys.h"
#include "ObitSwPower.h"
#include "ObitTableGC.h"
#include "ObitTableCD.h"
#include "ObitPrecess.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGainCal.c
 * ObitGainCal utility function definitions.
 */
/*---------------Private functions---------------------------*/
/* Read GC table and save */
static void GetGains (ObitUV *inData, olong GCVer, 
		      olong maxTerm, olong  numIF, olong numPoln,
		      olong  numAnt, 
		      ofloat *gains1, ofloat *sens1, 
		      ofloat *gains2, ofloat *sens2,
		      ObitErr *err);
/* Read CD table and save */
static void GetTCals (ObitUV *inData, olong CDVer, 
		      olong numIF, olong numPoln, olong  numAnt, 
		      ofloat *TCals1, ofloat *TCals2, ObitErr *err);
/* Calculate gains */
static void CalcGain (ObitTableRow *row, olong numIF, olong numPoln, olong maxTerm, 
		      gboolean doGain, gboolean doSens, gboolean doTsys, 
		      gboolean doOpac, gboolean doSwPwr,
		      ObitOpacity *Opac, ObitTsys *Tsys, ObitSwPower *SwPwr, 
		      odouble *freqs, ObitAntennaList **antList, ObitSource *source,
		      ofloat *TCals1, ofloat *Tcals2, 
		      ofloat *gains1, ofloat *sens1, ofloat *gains2, ofloat *sens2, 
		      ofloat *opac, 
		      ofloat *PwrDif1, ofloat *PwrSum1, ofloat *PwrGain1, 
		      ofloat *PwrDif2, ofloat *PwrSum2, ofloat *PwrGain2,
		      ofloat *gain1, ofloat *gain2, ObitErr *err);
/*---------------Public functions---------------------------*/
/**
 *  Calculate a CL or SN  table for inData
 *  A CL table is generated on inData for the times given in the iNdeX
 *  table with various optional amplitude calibrations.
 * \param inData UV data to be processed, info has following control parameters:
 * \li calInt   OBIT_float (1) Calibration interval in sec [def 30 sec].
 * \li calVer   OBIT_int   (1) CL table version, 0=>new [def 0]
 * \li doGain   OBIT_boo   (1) Apply antenna gain. [def FALSE]
 * \li doSens   OBIT_boo   (1) Apply nomiman sensitivity  [def FALSE]
 * \li GCVer    OBIT_int   (1) Gain table to use, 0=>highest [def 0]
 * \li doOpac   OBIT_boo   (1) Apply opacity correction? [def FALSE]
 * \li WXVer    OBIT_int   (1) Weather (WX) table to use, 0=>highest [def 0]
 * \li WXWeight OBIT_float (1) Weight of weather table in opacity [def 1]
 * \li doTsys   OBIT_boo   (1) Apply Tsys correction  [def FALSE]
 * \li TYVer    OBIT_int   (1) Tsys table to use, 0=>highest [def 0]
 * \li doSwPwr  OBIT_boo   (1) Apply EVLA switched power correction  [def FALSE]
 * \li SYVer    OBIT_int   (1) SY table to use, 0=>highest [def 0]
 * \li CDVer    OBIT_int   (1) CD table to use, 0=>highest [def 0]
 * \param doSN   If TRUE write SN table, else CL   
 * \param err    ObitError/message stack
 * \return Pointer to the newly created ObitTableCL object which is 
 *                 associated with outUV.  Table cast to ObitTable*
 */
ObitTable* ObitGainCalCalc (ObitUV *inData, gboolean doSN, ObitErr *err)
{
  ObitTable *outCal=NULL;
  ObitTableRow *row=NULL;
  gboolean doGain=FALSE, doSens=FALSE, doTsys=FALSE, doOpac=FALSE, doSwPwr=FALSE;
  olong calVer=0, GCVer=0, TYVer=0, WXVer=0, SYVer=0, CDVer=0;
  olong highVer, ver, iver;
  ofloat calInt=30.0;
  ObitTableAN    *ANTable=NULL;
  ObitTableSY    *SYTab=NULL;
  ObitTableTY    *TYTab=NULL;
  ObitTableSU    *SUTab=NULL;
  ObitTableNX    *NXTab=NULL;
  ObitTableNXRow *NXRow=NULL;
  ObitTableCL    *CLCal=NULL;
  ObitTableCLRow *CLRow=NULL;
  ObitTableSN    *SNCal=NULL;
  ObitTableSNRow *SNRow=NULL;
  ObitUVDesc     *desc=NULL;
  ObitSource     *source=NULL, *oneSource=NULL;
  ObitSourceList *souList=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat t0, endTime;
  ofloat lastTime=-1.0;
  olong iRow, oRow, i, ia, lrec;
  olong  numSubA, nfreq;
  oint numPoln=0, numIF=0, numTerm=0, numAnt=0, numOrb=0, numPCal=0;
  ObitAntennaList **antList=NULL;
  ObitOpacity *Opac=NULL;
  ObitTsys    *Tsys=NULL;
  ObitSwPower *SwPwr=NULL;
  ofloat      *gains1=NULL, *sens1=NULL, *gains2=NULL, *sens2=NULL;
  ofloat      *PwrDif1=NULL, *PwrSum1=NULL, *PwrGain1=NULL;
  ofloat      *PwrDif2=NULL, *PwrSum2=NULL, *PwrGain2=NULL;
  ofloat      *opac=NULL, *TCals1=NULL, *TCals2=NULL, *gain1=NULL, *gain2=NULL;
  odouble     *freqs=NULL;
  olong       maxTerm = 4;
  ObitIOCode  retCode;
  gchar *tname;
  gchar *routine = "ObitGainCalCalc";
 
   /* error checks */
  if (err->error) return outCal;
  g_assert (ObitUVIsA(inData));
  desc = inData->myDesc;

  /* Check that there is an NX table */
  highVer = ObitTableListGetHigh (inData->tableList, "AIPS NX");
  Obit_retval_if_fail((highVer>0), err, outCal,
		      "%s: NO iNdeX table found on %s", 
		      routine, inData->name);

  /* get control */
  ObitInfoListGetTest(inData->info, "calInt",  &type, dim, &calInt);
  calInt /= 86400.0;  /* to days */
  ObitInfoListGetTest(inData->info, "calVer",  &type, dim, &calVer);
  ObitInfoListGetTest(inData->info, "doGain",  &type, dim, &doGain);
  ObitInfoListGetTest(inData->info, "doSens",  &type, dim, &doSens);
  ObitInfoListGetTest(inData->info, "GCVer",   &type, dim, &GCVer);
  ObitInfoListGetTest(inData->info, "doOpac",  &type, dim, &doOpac);
  ObitInfoListGetTest(inData->info, "WXVer",   &type, dim, &WXVer);
  ObitInfoListGetTest(inData->info, "doSwPwr", &type, dim, &doSwPwr);
  ObitInfoListGetTest(inData->info, "SYVer",   &type, dim, &SYVer);
  ObitInfoListGetTest(inData->info, "CDVer",   &type, dim, &CDVer);

  /* open UV data to fully instantiate if not already open */
  if ((inData->myStatus==OBIT_Inactive) || (inData->myStatus==OBIT_Defined)) {
    retCode = ObitUVOpen (inData, OBIT_IO_ReadOnly, err);
    if (err->error) Obit_traceback_val (err, routine, inData->name, outCal);
  }
  lrec = inData->myDesc->lrec;
  t0 = -1.0e20;

  /* Init/Open index table */
  ver   = 1;
  tname = g_strconcat ("Index for: ", inData->name, NULL);
  NXTab = newObitTableNXValue(tname, (ObitData*)inData, &ver, OBIT_IO_ReadOnly,  
			      err);
  g_free (tname);
  ObitTableNXOpen (NXTab, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, outCal);
  /* Create Row */
  NXRow = newObitTableNXRow (NXTab);

  /* Create output */
  if (desc->jlocs>=0)  numPoln = MIN (2, desc->inaxes[desc->jlocs]);
  else                 numPoln = 1;
  if (desc->jlocif>=0) numIF = desc->inaxes[desc->jlocif];
  else                 numIF = 1;
  numTerm= 1;
  numAnt = inData->myDesc->numAnt[0];/* actually highest antenna number */
  tname  = g_strconcat ("Calibration for: ",inData->name, NULL);
  if (doSN) {
    SNCal = newObitTableSNValue(tname, (ObitData*)inData, &calVer, OBIT_IO_WriteOnly,  
				 numPoln, numIF, err);
    g_free (tname);
    outCal = (ObitTable*)SNCal;
    /* Create Row */
    SNRow = newObitTableSNRow (SNCal);
    row   = (ObitTableRow*)SNRow;
    
    /* Open table */
    ObitTableSNOpen (SNCal, OBIT_IO_WriteOnly, err) ;

    /* Attach row to output buffer */
    ObitTableSNSetRow (SNCal, SNRow, err);
    if (err->error) Obit_traceback_val (err, routine, inData->name, outCal);
    /* Set header values */
    SNCal->numAnt    = numAnt;  /* Max. antenna number */
  } else {
    CLCal = newObitTableCLValue(tname, (ObitData*)inData, &calVer, OBIT_IO_WriteOnly,  
				numPoln, numIF, numTerm, err);
    g_free (tname);
    outCal = (ObitTable*)CLCal;
    
    /* Create Row */
    CLRow = newObitTableCLRow (CLCal);
    row   = (ObitTableRow*)CLRow;
    
    /* Open table */
    ObitTableCLOpen (CLCal, OBIT_IO_WriteOnly, err);

    /* Attach row to output buffer */
    ObitTableCLSetRow (CLCal, CLRow, err);
    if (err->error) Obit_traceback_val (err, routine, inData->name, outCal);

    /* Set header values */
    CLCal->numAnt    = numAnt;  /* Max. antenna number */
  } /* end create CL table */


  /* Get additional information */
  /* Want Opacity? */
  if (doOpac) {
    opac  = g_malloc0(numIF*sizeof(ofloat));
    freqs = g_malloc0(numIF*sizeof(odouble));
    nfreq = numIF;
    for (i=0; i<numIF; i++) freqs[i] = desc->freqIF[i];
    Opac = ObitOpacityCreate ("Opacity", inData);
  }/* End init opacity */

  /* Want EVLA Switched power correction? */
  /* Check that tables available */  
  if (doSwPwr) {
    highVer = ObitTableListGetHigh (inData->tableList, "AIPS SY");
    if (highVer<=0) {
      doSwPwr = FALSE;
      Obit_log_error(err, OBIT_InfoWarn, 
		     "NO SY Table - cannot make EVLA switched power corrections");
    }
    highVer = ObitTableListGetHigh (inData->tableList, "AIPS CD");
    if (highVer<=0) {
      doSwPwr = FALSE;
      Obit_log_error(err, OBIT_InfoWarn, 
		     "NO CD Table - cannot make EVLA switched power corrections");
    }
  } /* end table testing */
  if (doSwPwr) {
    PwrDif1  = g_malloc0(numIF*sizeof(ofloat));
    PwrDif2  = g_malloc0(numIF*sizeof(ofloat));
    PwrSum1  = g_malloc0(numIF*sizeof(ofloat));
    PwrSum2  = g_malloc0(numIF*sizeof(ofloat));
    PwrGain1 = g_malloc0(numIF*sizeof(ofloat));
    PwrGain2 = g_malloc0(numIF*sizeof(ofloat));
    SYTab = newObitTableSYValue("SYTable", (ObitData*)inData, &SYVer, OBIT_IO_ReadOnly,  
				numIF, numPoln, err);
    SwPwr    = ObitSwPowerCreate ("SwPower", SYTab, inData, err);
    TCals1   = g_malloc0(numAnt*numIF*sizeof(ofloat));
    TCals2   = g_malloc0(numAnt*numIF*sizeof(ofloat));
    GetTCals(inData, CDVer, numIF, numPoln, numAnt, TCals1, TCals2, err);
    if (err->error) Obit_traceback_val (err, routine, inData->name, outCal);
  }/* End init Tsys */

  /* Want Tsys correction? */
  /* Check that table available */  
  if (doTsys) {
    highVer = ObitTableListGetHigh (inData->tableList, "AIPS TY");
    if (highVer<=0) {
      doTsys = FALSE;
      Obit_log_error(err, OBIT_InfoWarn, 
		     "NO TY Table - cannot make TSys corrections");
    }
  } /* end table check */
  if (doTsys) {
    TYTab = newObitTableTYValue("TYTable", (ObitData*)inData, &TYVer, OBIT_IO_ReadOnly,  
				numIF, numPoln, err);
    Tsys = ObitTsysCreate ("Tsys", TYTab, inData, err);
    if (err->error) Obit_traceback_val (err, routine, inData->name, outCal);
  }/* End init Tsys */

  /* Want Gain or sensitivity correction? */
  /* Check that table available */  
  if (doGain || doSens) {
    highVer = ObitTableListGetHigh (inData->tableList, "AIPS GC");
    if (highVer<=0) {
      doGain = doSens = FALSE;
      Obit_log_error(err, OBIT_InfoWarn, 
		     "NO GC Table - cannot make Gain or sensitivity corrections");
    }
  } /* end table check */
  if (doGain || doSens) {
    gains1 = g_malloc0(maxTerm*numIF*numAnt*sizeof(ofloat));
    sens1  = g_malloc0(numIF*numAnt*sizeof(ofloat));
    gains2 = g_malloc0(maxTerm*numIF*numAnt*sizeof(ofloat));
    sens2  = g_malloc0(numIF*numAnt*sizeof(ofloat));
    GetGains (inData, GCVer, maxTerm, numIF, numPoln, numAnt, 
	      gains1, sens1, gains2, sens2, err);
    if (err->error) Obit_traceback_val (err, routine, inData->name, outCal);
  }/* End init gain/sensitivity */
  ObitErrLog(err); /* show any  messages on err */

  /* Get antenna info */
  highVer = ObitTableListGetHigh (inData->tableList, "AIPS AN");
  numSubA = highVer;
  antList = g_malloc0(numSubA*sizeof(ObitAntennaList*));

  /* Loop over AN tables (subarrays) forming antenna lists*/
  for (iver=1; iver<=numSubA; iver++) {
    ver = iver;
    ANTable = newObitTableANValue ("AN table", (ObitData*)inData, 
				   &ver, OBIT_IO_ReadOnly, numIF, numOrb, numPCal, err);
    if (ANTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with AN table");
    antList[iver-1] = ObitTableANGetList (ANTable, err);
    if (err->error) Obit_traceback_val (err, routine, inData->name, outCal);
    
    ANTable = ObitTableANUnref(ANTable);   /* Cleanup */
  } /* End loop over subarrays */

  /* Source info */
  /* Is in->myData a single or multi source file? */
  highVer = ObitTableListGetHigh (inData->tableList, "AIPS SU");

  if (highVer>0) {
    /* Multisource, get Source List */
  /* Convert SU table into Source List */
    ver = 1;
    SUTab = newObitTableSUValue (inData->name, (ObitData*)inData, 
				   &ver, numIF, OBIT_IO_ReadOnly, err);
    if (SUTab) souList = ObitTableSUGetList (SUTab, err);
    if (err->error) Obit_traceback_val (err, routine, inData->name, outCal);
    SUTab = ObitTableSUUnref(SUTab);
  } else {
    /* single source - create oneSource */
    oneSource = newObitSource("Single");
    strncpy (oneSource->SourceName, inData->myDesc->object, 
	     MIN(20,UVLEN_VALUE));
    oneSource->equinox = inData->myDesc->equinox;
    oneSource->RAMean  = inData->myDesc->crval[inData->myDesc->jlocr];
    oneSource->DecMean = inData->myDesc->crval[inData->myDesc->jlocd];
    /* Compute apparent position */
    ObitPrecessUVJPrecessApp (inData->myDesc, oneSource);
  }

  /* gain work arrays */
  gain1 = g_malloc0(numIF*sizeof(ofloat));
  gain2 = g_malloc0(numIF*sizeof(ofloat));


 /* Initialize output row */
  if (doSN) {  /* SN Table */
    SNRow->Time   = 0.0;
    SNRow->TimeI  = 0.0;
    SNRow->SourID = 0;
    SNRow->antNo  = 0;
    SNRow->SubA   = 1;
    SNRow->FreqID = 0;
    SNRow->MBDelay1  = 0.0;
    /* IF dependent things */
    for (i=0; i<numIF; i++) {
      SNRow->Real1[i]   = 1.0;
      SNRow->Imag1[i]   = 0.0;
      SNRow->Rate1[i]   = 0.0;
      SNRow->Delay1[i]  = 0.0;
      SNRow->Weight1[i] = 1.0;
      SNRow->RefAnt1[i] = 0;
    }
    /* Multiple ppolarizations */
    if (numPoln>1) {
      SNRow->MBDelay1  = 0.0;
      /* IF dependent things */
      for (i=0; i<numIF; i++) {
	SNRow->Real2[i]   = 1.0;
	SNRow->Imag2[i]   = 0.0;
	SNRow->Rate2[i]   = 0.0;
	SNRow->Delay2[i]  = 0.0;
	SNRow->Weight2[i] = 1.0;
	SNRow->RefAnt2[i] = 0;
      }
    } /* end two poln */
  } else {  /* CL Table */
    CLRow->Time   = 0.0;
    CLRow->TimeI  = 0.0;
    CLRow->SourID = 0;
    CLRow->antNo  = 0;
    CLRow->SubA   = 1;
    CLRow->FreqID = 0;
    CLRow->GeoDelay[0]= 0.0;
    CLRow->IFR       = 0.0;
    CLRow->atmos     = 0.0;
    CLRow->Datmos    = 0.0;
    CLRow->MBDelay1  = 0.0;
    CLRow->clock1    = 0.0;
    CLRow->Dclock1   = 0.0;
    CLRow->dispers1  = 0.0;
    CLRow->Ddispers1 = 0.0;
    /* IF dependent things */
    for (i=0; i<numIF; i++) {
      CLRow->Real1[i]   = 1.0;
      CLRow->Imag1[i]   = 0.0;
      CLRow->Rate1[i]   = 0.0;
      CLRow->Delay1[i]  = 0.0;
      CLRow->Weight1[i] = 1.0;
      CLRow->RefAnt1[i] = 0;
    }
    /* Multiple ppolarizations */
    if (numPoln>1) {
      CLRow->clock2    = 0.0;
      CLRow->Dclock2   = 0.0;
      CLRow->dispers2  = 0.0;
      CLRow->Ddispers2 = 0.0;
      /* IF dependent things */
      for (i=0; i<numIF; i++) {
	CLRow->Real2[i]   = 1.0;
	CLRow->Imag2[i]   = 0.0;
	CLRow->Rate2[i]   = 0.0;
	CLRow->Delay2[i]  = 0.0;
	CLRow->Weight2[i] = 1.0;
	CLRow->RefAnt2[i] = 0;
      }
    } /* end two poln */
  } /* end init CLRow */

  /* loop looking at iNdeX table */
  for (iRow=1; iRow<=NXTab->myDesc->nrow; iRow++) {
    ObitTableNXReadRow (NXTab, iRow, NXRow, err);
    if (err->error) goto cleanup;
    if (NXRow->status==-1) continue;
    
    /* Multi or single source? */
    if (souList) source = souList->SUlist[NXRow->SourID-1];
    else source = oneSource;

    /* First time of scan */
    t0          = NXRow->Time - 0.5*NXRow->TimeI;
    lastTime    = t0;
    endTime     = NXRow->Time + 0.5*NXRow->TimeI;
    if (doSN) {
      SNRow->Time   = lastTime;
      SNRow->TimeI  = 0.0;
      SNRow->SourID = NXRow->SourID;
      SNRow->FreqID = NXRow->FreqID;
      SNRow->SubA   = MAX(1,NXRow->SubA);
    } else {
      CLRow->Time   = lastTime;
      CLRow->TimeI  = 0.0;
      CLRow->SourID = NXRow->SourID;
      CLRow->FreqID = NXRow->FreqID;
      CLRow->SubA   = MAX(1,NXRow->SubA);
    }
        
    /* Loop over antennas calculating gains, writing */
    for (ia=1; ia<=numAnt; ia++) {
      oRow = -1;
      if (doSN) SNRow->antNo = ia;
      else      CLRow->antNo = ia;

      /* Calculate gain */
      CalcGain (row, numIF, numPoln, maxTerm, doGain, doSens, doTsys, doOpac, doSwPwr,
		Opac, Tsys, SwPwr, freqs, antList, source,
		TCals1, TCals2, gains1, sens1, gains2, sens2, 
		opac, PwrDif1, PwrSum1, PwrGain1, PwrDif2, PwrSum2, PwrGain2,
		gain1, gain2, err);
      /* Write */
      if (doSN) ObitTableSNWriteRow (SNCal, oRow, SNRow, err);
      else      ObitTableCLWriteRow (CLCal, oRow, CLRow, err);
      if (err->error) goto cleanup;
    } /* end initial antenna loop */

   /* Loop over scan */
    lastTime += calInt;
    while (lastTime<=endTime) {
      if (doSN) SNRow->Time   = lastTime;
      else      CLRow->Time   = lastTime;
      
      /* Loop over antennas calculating gains, writing */
      for (ia=1; ia<=numAnt; ia++) {
	oRow = -1;
	if (doSN) SNRow->antNo = ia;
	else      CLRow->antNo = ia;
	
	/* Calculate gain */
	CalcGain (row, numIF, numPoln, maxTerm, doGain, doSens, doTsys, doOpac, doSwPwr,
		  Opac, Tsys, SwPwr, freqs, antList, source,
		  TCals1, TCals2, gains1, sens1, gains2, sens2, 
		  opac, PwrDif1, PwrSum1, PwrGain1, PwrDif2, PwrSum2, PwrGain2,
		  gain1, gain2, err);
	/* Write */
	if (doSN) ObitTableSNWriteRow (SNCal, oRow, SNRow, err);
	else      ObitTableCLWriteRow (CLCal, oRow, CLRow, err);
	if (err->error) goto cleanup;
      } /* end  antenna loop */
      lastTime += calInt;
    } /* end loop over scan */

    /* End of scan */
    if (doSN) SNRow->Time   = endTime;
    else      CLRow->Time   = endTime;

    /* Loop over antennas calculating gains, writing */
    for (ia=1; ia<=numAnt; ia++) {
      oRow = -1;
      if (doSN) SNRow->antNo = ia;
      else      CLRow->antNo = ia;
      
      /* Calculate gain */
      CalcGain (row, numIF, numPoln, maxTerm, doGain, doSens, doTsys, doOpac, doSwPwr,
		Opac, Tsys, SwPwr, freqs,  antList, source,
		TCals1, TCals2, gains1, sens1, gains2, sens2, 
		opac, PwrDif1, PwrSum1, PwrGain1, PwrDif2, PwrSum2, PwrGain2,
		gain1, gain2, err);
      /* Write */
      if (doSN) ObitTableSNWriteRow (SNCal, oRow, SNRow, err);
      else      ObitTableCLWriteRow (CLCal, oRow, CLRow, err);
      if (err->error) goto cleanup;
    } /* end final antenna loop */
  } /* End scan loop */

  /* Close UV data */
  ObitUVClose(inData, err);
  if (err->error) goto cleanup;

  /* Close cal table */
  if (doSN) ObitTableSNClose (SNCal, err);
  else      ObitTableCLClose (CLCal, err);
  if (err->error) goto cleanup;
  
  /* Close NX table */
  retCode = ObitTableNXClose (NXTab, err);
  if (err->error) goto cleanup;

  /* Cleanup */
 cleanup:
  Opac   = ObitOpacityUnref(Opac);
  Tsys   = ObitTsysUnref(Tsys);
  SwPwr  = ObitSwPowerUnref(SwPwr);
  SYTab  = ObitTableSYUnref(SYTab);
  TYTab  = ObitTableTYUnref(TYTab);
  NXTab  = ObitTableNXUnref(NXTab);
  NXRow  = ObitTableNXRowUnref(NXRow);
  CLRow  = ObitTableCLRowUnref(CLRow);
  SNRow  = ObitTableCLRowUnref(SNRow);
  oneSource = ObitSourceUnref(oneSource);
  souList   = ObitSourceListUnref(souList);
  if (gain1)    g_free(gain1);
  if (gain2)    g_free(gain2);
  if (gains1)   g_free(gains1);
  if (gains2)   g_free(gains2);
  if (sens1)    g_free(sens1);
  if (sens2)    g_free(sens2);
  if (opac)     g_free(opac);
  if (freqs)    g_free(freqs);
  if (TCals1)   g_free(TCals1);
  if (TCals2)   g_free(TCals2);
  if (PwrDif1)  g_free(PwrDif1);
  if (PwrDif2)  g_free(PwrDif2);
  if (PwrSum1)  g_free(PwrSum1);
  if (PwrSum2)  g_free(PwrSum2);
  if (PwrGain1) g_free(PwrGain1);
  if (PwrGain2) g_free(PwrGain2);
  if (antList) {
    for (i=0; i<numSubA; i++) 
      antList[i] = ObitAntennaListUnref(antList[i]);
    g_free(antList);
  }
  if (err->error) Obit_traceback_val (err, routine, inData->name, outCal);

  return outCal;
} /* end ObitGainCalCalc */

/*---------------Private functions---------------------------*/
/**
 *  Read GC table into an array
 *  Only deals with type=2, polynomial in ZA
 * \param inData  UV data to be processed, info has following control parameters:
 * \param GCVer   Gain table version to use, 0=>highest 
 * \param maxTerm Maximum number of polynomial terms
 * \param numIF   Number of IFs
 * \param numPoln Number of polarizations
 * \param numAnt  Maximum 1-rel antenna number
 * \param gains1  [out] gain polynomial coefs, poln 1
 *                Unused terms fblanked,, all terms blanked -> flat curve
 *                From fastest to slowest, term, IF, Ant
 * \param sens1   [out] Nominal sensitivity (K/Jy) per antenna, poln 1
 *                From fastest to slowest, IF, Ant
 * \param gains2  [out] gain polynomial coefs, poln 2
 *                Unused terms fblanked,
 *                From fastest to slowest, term, IF, Ant
 * \param sens2   [out] Nominal sensitivity (K/Jy) per antenna, poln 2
 *                From fastest to slowest, IF, Ant
 * \param err     ObitError/message stack
 */
static void GetGains (ObitUV *inData, olong GCVer, 
		      olong maxTerm, olong  numIF, olong numPoln,
		      olong  numAnt, 
		      ofloat *gains1, ofloat *sens1,
		      ofloat *gains2, ofloat *sens2,
		      ObitErr *err)
{
  ObitTableGC    *GCTab=NULL;
  ObitTableGCRow *GCRow=NULL;
  olong i, iIF, n, iRow, nvalPAnt, nvalPIF, indx;
  oint numBand=0, numPol=0, numTabs=0;
  ofloat fblank = ObitMagicF();
  gchar *routine = "ObitGainCal:GetGains";

  /* Previous error? */
  if(err->error) return;

  /* Initialize output */
  n = maxTerm*numAnt*numIF;
  for (i=0; i<n; i++) gains1[i] = fblank;
  for (i=0; i<n; i++) gains2[i] = fblank;
  n =  numAnt*numIF;
  for (i=0; i<n; i++) sens1[i]  = 0.0;
  for (i=0; i<n; i++) sens2[i]  = 0.0;

  /* How many of various things? for gains */
  nvalPAnt  = maxTerm*numIF;         /* Values per antenna */
  nvalPIF   = maxTerm;               /* Values per IF */

  /* Make table, open */
  GCTab = 
    newObitTableGCValue (inData->name, (ObitData*)inData, &GCVer, OBIT_IO_ReadOnly, 
			 numBand, numPol, numTabs, err);
  ObitTableGCOpen (GCTab, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Get actual max number of gain entries */
  numTabs = GCTab->numTabs;
  
  GCRow  = newObitTableGCRow (GCTab);  /* Make row */

  /* Loop over table */
  for (iRow=1; iRow<=GCTab->myDesc->nrow; iRow++) {
    ObitTableGCReadRow (GCTab, iRow, GCRow, err);
    if (err->error) goto cleanup;
    if (GCRow->status==-1) continue;
    /* Save gain */
    /* Loop over IF */
    for (iIF=0; iIF<numIF; iIF++) {
      indx = (GCRow->antennaNo-1)*nvalPAnt + iIF*nvalPIF;
      for (i=0; i<GCRow->NTerm1[iIF]; i++) 
	gains1[indx+i] = GCRow->gain1[iIF*numTabs+i];
      /* Second Poln */
      if (numPoln>1) {
        for (i=0; i<GCRow->NTerm2[iIF]; i++) 
	  gains2[indx+i] = GCRow->gain2[iIF*numTabs+i];
      } /* end 2nd poln */
    } /* end loop over IF */

    /* Save sens */
    /* Loop over IF */
    for (iIF=0; iIF<numIF; iIF++) {
      indx = (GCRow->antennaNo-1)*numIF + iIF;
      sens1[indx] = GCRow->sens1[iIF];
      /* Second Poln */
      if (numPoln>1) {
	sens2[indx] = GCRow->sens2[iIF];
      } /* end sens second poln */
    } /* end loop over IF */
  } /* end loop over table */
  
  /* Check for no gain curve */
  for (i=0; i<numAnt; i++) {
    for (iIF=0; iIF<numIF; iIF++) {
      indx = i*nvalPAnt + iIF*nvalPIF;
      if (gains1[indx]==fblank) gains1[indx] = 1.0;  /* Flatliner */
      if ((numPoln>1) && (gains2[indx]==fblank)) gains2[indx] = 1.0;  /* Flatliner */
    }
  }

  
  /* Cleanup */
 cleanup:
  ObitTableGCClose (GCTab, err);
  GCRow = ObitTableGCRowUnref(GCRow);
  GCTab = ObitTableGCUnref(GCTab);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
} /* end GetGains */

/**
 *  Read CD (EVLA Tcal data) table into an array
 * \param inData  UV data to be processed, info has following control parameters:
 * \param CDVer   Table version to use, 0=>highest 
 * \param numIF   Number of IFs
 * \param numPoln Number of polarizations
 * \param numAnt  Maximum 1-rel antenna number
 * \param TCals1 [out] Tcal poln 1, per IF, per antenna (K)
 * \param TCals2 [out] Tcal poln 2, per IF, per antenna (K)
 *                fblanked if not available
 * \param err     ObitError/message stack
 */
static void GetTCals (ObitUV *inData, olong CDVer, 
		      olong numIF, olong numPoln, olong  numAnt, 
		      ofloat *TCals1, ofloat *TCals2, ObitErr *err)
{
  ObitTableCD    *CDTab=NULL;
  ObitTableCDRow *CDRow=NULL;
  olong i, iIF, n, iRow, indx;
  oint nIF=(oint)numIF, nPol=(oint)numPoln;
  gchar *routine = "ObitGainCal:GetTCals";

  /* Previous error? */
  if(err->error) return;

  /* Initialize output */
  n = numAnt*numIF;
  for (i=0; i<n; i++) TCals1[i] = -1.0;
  for (i=0; i<n; i++) TCals2[i] = -1.0;

  /* Make table, open */
  CDTab = 
    newObitTableCDValue (inData->name, (ObitData*)inData, &CDVer, OBIT_IO_ReadOnly, 
			 nIF, nPol, err);
  ObitTableCDOpen (CDTab, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  CDRow  = newObitTableCDRow (CDTab);  /* Make row */

  /* Loop over table */
  for (iRow=1; iRow<=CDTab->myDesc->nrow; iRow++) {
    ObitTableCDReadRow (CDTab, iRow, CDRow, err);
    if (err->error) goto cleanup;
    if (CDRow->status==-1) continue;
    /* Save  */
    /* Loop over IF */
    for (iIF=0; iIF<numIF; iIF++) {
      indx = (CDRow->antennaNo-1)*numIF + iIF;
      TCals1[indx] = CDRow->TCal1[iIF];
      /* Second Poln */
      if (numPoln>1) {
	TCals2[indx] = CDRow->TCal2[iIF];
      } /* end 2nd poln */
    } /* end loop over IF */

  } /* end loop over table */
  
  /* Cleanup */
 cleanup:
  ObitTableCDClose (CDTab, err);
  CDRow = ObitTableCDRowUnref(CDRow);
  CDTab = ObitTableCDUnref(CDTab);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
} /* end GetTCals */

/**
 *  Calculate Gain and apply to CL table row
 * \param row      CL or SN table row to update
 * \param numIF    Number of IFs
 * \param numPoln  Number of polarizations
 * \param maxTerm  Number of antenna gain terms
 * \param doGain   Apply antenna gain?
 * \param doSens   Apply nomimal sensitivity?
 * \param doOpac   Apply opacity correction?
 * \param doTsys   Apply Tsys correction?
 * \param doSwPwr  Apply EVLA switched power correction?
 * \param Opac     Opacity calculator
 * \param Tsys     TSys interpolator
 * \param SwPwr    Switched power interpolator
 * \param freqs    IF freqyency array (Hz)
 * \param antList  Array of antenna lists.
 * \param source   Source of interest
 * \param TCals1   Tcal poln 1, per IF, per antenna (K)
 * \param TCals2   Tcal poln 2, per IF, per antenna (K
 * \param gains1   gain polynomial coefs, poln 1
 * \param sens1    Nominal sensitivity (K/Jy) per antenna, poln 1
 * \param gains2   gain polynomial coefs, poln 2
 * \param sens2    Nominal sensitivity (K/Jy) per antenna, poln 2
 * \param opac     Work array
 * \param PwrDif1  Work array
 * \param PwrSum1  Work array
 * \param PwrGain1 Work array
 * \param PwrDif2  Work array
 * \param PwrSum2  Work array
 * \param PwrGain2 Work array
 * \param gain1    Work array
 * \param gain2    Work array
 * \param err      ObitError/message stack
 */
static void CalcGain (ObitTableRow *row, olong numIF, olong numPoln, olong maxTerm, 
		      gboolean doGain, gboolean doSens, gboolean doTsys, 
		      gboolean doOpac, gboolean doSwPwr,
		      ObitOpacity *Opac, ObitTsys *Tsys, ObitSwPower *SwPwr, 
		      odouble *freqs, ObitAntennaList **antList, ObitSource *source,
		      ofloat *TCals1, ofloat *TCals2, 
		      ofloat *gains1, ofloat *sens1, ofloat *gains2, ofloat *sens2, 
		      ofloat *opac, 
		      ofloat *PwrDif1, ofloat *PwrSum1, ofloat *PwrGain1, 
		      ofloat *PwrDif2, ofloat *PwrSum2, ofloat *PwrGain2,
		      ofloat *gain1, ofloat *gain2,ObitErr *err )
{
  ObitTableSNRow *SNRow=NULL;
  ObitTableCLRow *CLRow=NULL;
  gboolean isSN;
  ofloat sum, za, elev, g, arg;
  ofloat fblank = ObitMagicF();
  ofloat *Tsys1=PwrSum1, *Tsys2=PwrSum2;   /* Alias work array */
  olong i, iIF, indx, Ant, SubA, SourID, FreqID;
  ofloat Time;
  gchar *routine = "ObitGainCal:CalcGain";

  /* initial values */
  for (i=0; i<numIF; i++) gain1[i] = gain2[i] = 1.0;

  /* Which type of table */
  isSN = ObitTableSNIsA (row);
  if (isSN) {
    SNRow  = (ObitTableSNRow*)row;
    Time   = SNRow->Time;
    Ant    = SNRow->antNo;
    SubA   = MAX (1, SNRow->SubA);
    SourID = SNRow->SourID;
    FreqID = SNRow->FreqID;
  } else {
    CLRow = (ObitTableCLRow*)row;
    Time   = CLRow->Time;
    Ant    = CLRow->antNo;
    SubA   = MAX (1, CLRow->SubA);
    SourID = CLRow->SourID;
    FreqID = CLRow->FreqID;
  }


  /* Antenna gain */
  if (doGain) {
    /* Get elevation */
    elev = ObitAntennaListElev(antList[SubA-1], Ant, Time, source);
    za = 90.0 - elev*RAD2DG;   /* Zenith angle */

    /* Loop over IFs */
    for (iIF=0; iIF<numIF; iIF++) {
      indx = (Ant-1)*maxTerm*numIF + iIF*maxTerm;
      /* First poln */
      sum = gains1[indx];
      arg = za;
      for (i=1; i<maxTerm; i++) {
	if (gains1[indx+i]==fblank) break;
	sum += arg*gains1[indx+i];
	arg  *= za;
      }
      gain1[iIF] /= sum;   /* correction inverse of gain */
      /* Second poln */
      if (numPoln>1) {
	sum = gains2[indx];
	arg = za;
	for (i=1; i<maxTerm; i++) {
	  if (gains2[indx+i]==fblank) break;
	  sum += arg*gains2[indx+i];
	  arg  *= za;
	}
	gain2[iIF] /= sum;  /* correction inverse of gain */
      }
    } /* end IF loop */
  }/* end antenna gain */

  /* Nomimal sensitivity */
  if (doSens) {
    indx = (Ant-1)*numIF;
    for (iIF=0; iIF<numIF; iIF++) {
      gain1[iIF] *= sens1[indx+iIF];
      if (numPoln>1) gain2[iIF] *= sens2[indx+iIF];
    }
  } /* end  Nominal sensitivity */

  /* Tsys */
  if (doTsys) {
    ObitTsysReport (Tsys, Time, Ant, SubA, Tsys1, Tsys2, err);
    if (err->error) Obit_traceback_msg (err, routine, Tsys->name);
    for (iIF=0; iIF<numIF; iIF++) {
      gain1[iIF] *= sqrt(Tsys1[iIF]);
      if (numPoln>1) gain2[iIF] *= sqrt(Tsys1[iIF]);
    }
  } /* end Tsys */

  /* Opacity */
  if (doOpac) {
    ObitOpacityCalc (Opac, Time, numIF, freqs, Ant, SubA, SourID, opac, err);
    if (err->error) Obit_traceback_msg (err, routine, Opac->name);
    for (iIF=0; iIF<numIF; iIF++) {
      g = opac[iIF];
      gain1[iIF] *= g;
      if (numPoln>1) gain2[iIF] *= g;
    }
  } /* end Opacity */

  /* EVLA switched power as per R. Perley EVLA memo 145 */
  if (doSwPwr) {
    ObitSwPowerReport (SwPwr, Time, Ant, SubA, PwrDif1, PwrSum1, PwrGain1,
		       PwrDif2, PwrSum2, PwrGain2, err);
    if (err->error) Obit_traceback_msg (err, routine, SwPwr->name);
    indx = (Ant-1)*numIF;
    for (iIF=0; iIF<numIF; iIF++) {
      if (fabs(PwrDif1[iIF])>0.0)
	gain1[iIF] *= sqrt (TCals1[indx+iIF] / fabs(PwrDif1[iIF]));
      if (numPoln>1) {
	if (fabs(PwrDif2[iIF])>0.0)
	gain2[iIF] *= sqrt (TCals2[indx+iIF] / fabs(PwrDif2[iIF]));
      }
    } /* end IF loop */
  } /* end EVLA switched power */
  
  /* apply */
  if (isSN) { /* SN table */
    for (iIF=0; iIF<numIF; iIF++) {
      SNRow->Real1[iIF] = gain1[iIF];
      if (numPoln>1) SNRow->Real2[iIF] = gain2[iIF];
    }
  } else {   /* CL table */
    for (iIF=0; iIF<numIF; iIF++) {
      CLRow->Real1[iIF] = gain1[iIF];
      if (numPoln>1) CLRow->Real2[iIF] = gain2[iIF];
    }
  }
} /* end CalcGain  */
