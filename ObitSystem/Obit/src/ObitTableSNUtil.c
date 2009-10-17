/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005,2009                                          */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include <math.h>
#include "ObitTableSNUtil.h"
#include "ObitTableANUtil.h"
#include "ObitTableSUUtil.h"
#include "ObitUVUtil.h"
#include "ObitPrecess.h"
#define MAXANT    300    /* Maximum number of antennas */

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableSNUtil.c
 * ObitTableSN class utility function definitions.
 */

/*---------------Private function prototypes----------------*/
/** Private: Sinc function. */
static ofloat Sinc (ofloat arg);

/** Private: Write SN table for ZeroFR */
static void SetRow (ObitAntennaList *AntList, ObitSourceList *SouList, 
		    odouble ArrLong, ObitUVDesc *desc,
		    odouble *RAR, odouble *DecR, ofloat delTime,
		    olong maxant, gboolean  gotAnt[MAXANT], gboolean invert,
		    ObitTableSNRow *row, ObitTableSN *outCal, ObitErr *err);

/*----------------------Public functions---------------------------*/

/**
 * Copies SN tables from inUV to outUV with selection in inUV
 * If calibration is selected on inUV, no tables are copied
 * \param inUV     Input UV to copy from
 * \param outUV    Output UV to copy to
 * \param err      ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableSNSelect (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableSN    *inTab=NULL, *outTab=NULL;
  ObitTableSNRow *inRow=NULL, *outRow=NULL;
  ObitInfoType type;
  olong itemp, iif, oif;
  olong highSNver, iSNver, inSNRow, outSNRow;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  union ObitInfoListEquiv InfoReal; 
  oint numPol, numIF;
  gboolean wanted;
  gchar *SNType = "AIPS SN";
  gchar *routine = "ObitTableSNSelect";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));

  /* Calibration selected? */
  type = OBIT_float; InfoReal.flt = 0.0;
  ObitInfoListGetTest(inUV->info, "doCalib", &type, (gint32*)dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  if (itemp>0) return  OBIT_IO_OK;  /* Yep, don't copy */

  /* Fully instantiate UV files */
  ObitUVFullInstantiate (inUV, TRUE, err);
  if (err->error )Obit_traceback_val (err, routine, inUV->name, retCode);
  ObitUVFullInstantiate (outUV, FALSE, err);
  if (err->error )Obit_traceback_val (err, routine, outUV->name, retCode);

  /* How many SN tables  */
  highSNver = ObitTableListGetHigh (inUV->tableList, SNType);

  /* Are there any? */
  if (highSNver <= 0) return OBIT_IO_OK;

  /* Loop over SN tables */
  for (iSNver=1; iSNver<=highSNver; iSNver++) {

    /* Get input table */
    numPol = 0;
    numIF = 0;
    inTab = 
      newObitTableSNValue (inUV->name, (ObitData*)inUV, &iSNver, OBIT_IO_ReadOnly, 
			   numPol, numIF, err);
    if (err->error) Obit_traceback_val (err, routine, inTab->name, retCode);
    /* Find it */
    if (inTab==NULL) continue;  /* No keep looping */

    /* Open input table */
    retCode = ObitTableSNOpen (inTab, OBIT_IO_ReadOnly, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inTab->name, retCode);

    /* Delete any old output table */
    retCode = ObitDataZapTable ((ObitData*)outUV, SNType, iSNver, err);
    if (err->error) Obit_traceback_val (err, routine, outUV->name, retCode);
 
    /* Create output table */
    numPol = MIN (2, inUV->mySel->numberPoln);
    numIF  = inUV->mySel->numberIF;
    outTab = 
      newObitTableSNValue (outUV->name, (ObitData*)outUV, &iSNver, OBIT_IO_WriteOnly, 
			   numPol, numIF, err);
    if (err->error) Obit_traceback_val (err, routine, outUV->name, retCode);
    /* Create it? */
    Obit_retval_if_fail((outTab!=NULL), err, retCode,
			"%s: Could not create SN table %d for %s", 
			routine, iSNver, outTab->name);

  /* Open output table */
  retCode = ObitTableSNOpen (outTab, OBIT_IO_WriteOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inTab->name, retCode);
  
    /* Update header info */
    outTab->revision  = inTab->revision;
    outTab->numAnt    = inTab->numAnt;
    outTab->numNodes  = inTab->numNodes;
    outTab->mGMod     = inTab->mGMod;;
    outTab->isApplied = inTab->isApplied;

    /* Set rows */
    inRow  = newObitTableSNRow (inTab);
    outRow = newObitTableSNRow (outTab);
    ObitTableSNSetRow (outTab, outRow, err);
    if (err->error) Obit_traceback_val (err, routine, outTab->name, retCode);

    /* Loop over table copying selected data */
    outSNRow = -1;
    for (inSNRow=1; inSNRow<=inTab->myDesc->nrow; inSNRow++) {
      retCode = ObitTableSNReadRow (inTab, inSNRow, inRow, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, inUV->name, retCode);
      if (inRow->status==-1) continue;
  
      /* Want this one? */
      wanted = ((inRow->Time >= inUV->mySel->timeRange[0]) && 
		(inRow->Time <= inUV->mySel->timeRange[1]));
      wanted = wanted && ((inUV->mySel->FreqID <= 0) ||
			  (inRow->FreqID==inUV->mySel->FreqID));
      wanted = wanted && ((inUV->mySel->SubA <= 0) ||
			  (inRow->SubA==inUV->mySel->SubA));
      wanted = wanted && ObitUVSelWantAnt(inUV->mySel, inRow->antNo);
      wanted = wanted && ObitUVSelWantSour(inUV->mySel, inRow->SourID);
      if (!wanted) continue;

      /* Copy selected data */
      outRow->Time     = inRow->Time;
      outRow->TimeI    = inRow->TimeI;
      outRow->SourID   = inRow->SourID;
      outRow->antNo    = inRow->antNo;
      outRow->SubA     = inRow->SubA;
      outRow->FreqID   = inRow->FreqID;
      outRow->IFR      = inRow->IFR;
      outRow->NodeNo   = inRow->NodeNo;
      outRow->MBDelay1 = inRow->MBDelay1;
      oif = 0;
      for (iif=inUV->mySel->startIF-1; 
	   iif<inUV->mySel->startIF+inUV->mySel->numberIF-1;
	   iif++) {
	     outRow->Real1[oif]   = inRow->Real1[iif];
	     outRow->Imag1[oif]   = inRow->Imag1[iif];
	     outRow->Delay1[oif]  = inRow->Delay1[iif];
	     outRow->Rate1[oif]   = inRow->Rate1[iif];
	     outRow->Weight1[oif] = inRow->Weight1[iif];
	     outRow->RefAnt1[oif] = inRow->RefAnt1[iif];
	     oif++;
	   }
      if (numPol>1) {
	outRow->MBDelay2 = inRow->MBDelay2;
	oif = 0;
	for (iif=inUV->mySel->startIF-1; 
	     iif<inUV->mySel->startIF+inUV->mySel->numberIF-1;
	     iif++) {
	  outRow->Real2[oif]   = inRow->Real2[iif];
	  outRow->Imag2[oif]   = inRow->Imag2[iif];
	  outRow->Delay2[oif]  = inRow->Delay2[iif];
	  outRow->Rate2[oif]   = inRow->Rate2[iif];
	  outRow->Weight2[oif] = inRow->Weight2[iif];
	  outRow->RefAnt2[oif] = inRow->RefAnt2[iif];
	  oif++;
	}
      } /* End second poln */

      retCode = ObitTableSNWriteRow (outTab, outSNRow, outRow, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, inUV->name, retCode);
    } /* end loop over rows */
    
    /* Close tables */
    retCode = ObitTableSNClose (inTab, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inTab->name, retCode);
    retCode = ObitTableSNClose (outTab, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, outTab->name, retCode);
 
    /* release table objects */
    inTab  = ObitTableSNUnref(inTab);
    outTab = ObitTableSNUnref(outTab);

    /* release row objects */
    inRow  = ObitTableSNRowUnref(inRow);
    outRow = ObitTableSNRowUnref(outRow);
  } /* end loop over tables */

  return retCode;
} /* end ObitTableSNSelect */

/**
 * Invert gains in an AIPS SolutioN table.
 * This procedure is used in reverting calibration in a data set.
 * The amplitudes are the reciprocal of the input and the phases, 
 * delays, rates  are the negative of the input.
 * \param inSN     The Table to invert
 * \param outData  The data object to which to attach the output table
 * \param outVer   Output SN table version, 0=> create new
 * \param doRepl   If True replace failed solutions with (1,0)
 * \param err      ObitErr error stack.
 * \return pointer to the new object.
 */
ObitTableSN* ObitTableSNUtilInvert (ObitTableSN *inSN, ObitData *outData, olong *outVer, 
				    gboolean doRepl, ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  ObitTableSN    *outSN=NULL;
  ObitTableSNRow *inRow=NULL, *outRow=NULL;
  ofloat amp, phase, fblank = ObitMagicF();
  ofloat replReal, replImag, replDelay, replRate, replWeight;
  olong inSNRow, outSNRow, oif, iif;
  gchar *routine = "ObitTableSNUtilInvert";
  
  /* error checks */
  if (err->error) return outSN;
  g_assert (ObitTableSNIsA(inSN));
  g_assert (ObitDataIsA(outData));

  /* Replacement for flagged solutions */
  if (doRepl) {
    replReal   = 1.0;
    replImag   = 0.0;
    replDelay  = 0.0;
    replRate   = 0.0;
    replWeight = 1.0;
  } else { /* leave flagged */
    replReal   = fblank;
    replImag   = fblank;
    replDelay  = fblank;
    replRate   = fblank;
    replWeight = 0.0;
  }
  
  /* Open input table */
  iretCode = ObitTableSNOpen (inSN, OBIT_IO_ReadOnly, err);
  if ((iretCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inSN->name, outSN);
  
  /* Create output table */
  outSN = 
    newObitTableSNValue (outData->name, outData, outVer, OBIT_IO_WriteOnly, 
			 inSN->numPol, inSN->numIF, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, outSN);
  /* Create it? */
  Obit_retval_if_fail((outSN!=NULL), err, outSN,
		      "%s: Could not create SN table %d for %s", 
		      routine, *outVer, outData->name);
  
  /* Open output table */
  oretCode = ObitTableSNOpen (outSN, OBIT_IO_WriteOnly, err);
  if ((oretCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inSN->name, outSN);
  
  /* Update header info */
  outSN->revision  = inSN->revision;
  outSN->numAnt    = inSN->numAnt;
  outSN->numNodes  = inSN->numNodes;
  if (inSN->mGMod>0.0) outSN->mGMod = 1.0 / inSN->mGMod;
  else outSN->mGMod = 1.0;
  outSN->isApplied  = FALSE;
  
  /* Set rows */
  inRow  = newObitTableSNRow (inSN);
  outRow = newObitTableSNRow (outSN);
  ObitTableSNSetRow (outSN, outRow, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, outSN);

    /* Loop over table inverting data */
    outSNRow = -1;
    for (inSNRow=1; inSNRow<=inSN->myDesc->nrow; inSNRow++) {
      iretCode = ObitTableSNReadRow (inSN, inSNRow, inRow, err);
      if ((iretCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, inSN->name, outSN);
      if (inRow->status==-1) continue;
  
      /* Copy/invert data */
      outRow->Time     = inRow->Time;
      outRow->TimeI    = inRow->TimeI;
      outRow->SourID   = inRow->SourID;
      outRow->antNo    = inRow->antNo;
      outRow->SubA     = inRow->SubA;
      outRow->FreqID   = inRow->FreqID;
      outRow->IFR      = inRow->IFR;
      outRow->NodeNo   = inRow->NodeNo;
      outRow->MBDelay1 = -inRow->MBDelay1;
      oif = 0;
      for (iif=0; iif<inSN->numIF; iif++) {
	if ((inRow->Weight1[iif]>0.0) && 
	    (inRow->Real1[iif]!=fblank)  && (inRow->Imag1[iif]!=fblank)) {
	  amp = sqrt (inRow->Real1[iif]*inRow->Real1[iif] + 
		      inRow->Imag1[iif]*inRow->Imag1[iif]);
	  if (amp>0.0) {
	    amp = 1.0 / amp;
	    phase = -atan2 (inRow->Imag1[iif], inRow->Real1[iif]);
	  } else {
	    phase = 0.0;
	  }
	  outRow->Real1[oif]   = amp * cos (phase);
	  outRow->Imag1[oif]   = amp * sin (phase);
	  outRow->Delay1[oif]  = -inRow->Delay1[iif];
	  outRow->Rate1[oif]   = -inRow->Rate1[iif];
	  outRow->Weight1[oif] =  inRow->Weight1[iif];
	} else { /* blanked */
	  outRow->Real1[oif]   = replReal;
	  outRow->Imag1[oif]   = replImag;
	  outRow->Delay1[oif]  = replDelay;
	  outRow->Rate1[oif]   = replRate;
	  outRow->Weight1[oif] = replWeight;
	}
	outRow->RefAnt1[oif] = inRow->RefAnt1[iif];
	oif++;
      }
      if (inSN->numPol>1) {
	outRow->MBDelay2 = -inRow->MBDelay2;
	oif = 0;
	for (iif=0; iif<inSN->numIF; iif++) {
	  if ((inRow->Weight2[iif]>0.0) && 
	      (inRow->Real2[iif]!=fblank) && (inRow->Imag2[iif]!=fblank)) {
	    amp = sqrt (inRow->Real2[iif]*inRow->Real2[iif] + 
			inRow->Imag2[iif]*inRow->Imag2[iif]);
	    if (amp>0.0) {
	      amp = 1.0 / amp;
	      phase = -atan2 (inRow->Imag2[iif], inRow->Real2[iif]);
	    } else {
	      phase = 0.0;
	    }
	    outRow->Real2[oif]   =  amp * cos (phase);
	    outRow->Imag2[oif]   =  amp * sin (phase);
	    outRow->Delay2[oif]  = -inRow->Delay2[iif];
	    outRow->Rate2[oif]   = -inRow->Rate2[iif];
	    outRow->Weight2[oif] =  inRow->Weight2[iif];	
	  } else { /* blanked */
	    outRow->Real2[oif]   = replReal;
	    outRow->Imag2[oif]   = replImag;
	    outRow->Delay2[oif]  = replDelay;
	    outRow->Rate2[oif]   = replRate;
	    outRow->Weight2[oif] = replWeight;
	  }
	  outRow->RefAnt2[oif] = inRow->RefAnt2[iif];
	  oif++;
	}
      } /* End second poln */

      oretCode = ObitTableSNWriteRow (outSN, outSNRow, outRow, err);
      if ((oretCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, outSN->name, outSN);
    } /* end loop over rows */
    
    /* Close tables */
    iretCode = ObitTableSNClose (inSN, err);
    if ((iretCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inSN->name, outSN);
    oretCode = ObitTableSNClose (outSN, err);
    if ((oretCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, outSN->name, outSN);
 
    /* release row objects */
    inRow  = ObitTableSNRowUnref(inRow);
    outRow = ObitTableSNRowUnref(outRow);

    /* Done */
    return outSN;
} /* end ObitTableSNUtilInvert */

#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif
/**
 * Create SN table that will counter-rotate the data to a zero fringe rate
 * After this operation, all terrestial sources should be constant.
 * Weights reflect the effect of the difference in fringe rate and delay
 * for the integration time and observed bandwidth.
 * \param inUV     Input UV data. Control parameters:
 * \li "solInt"    OBIT_float (1,1,1) Entry interval in days [def 10 sec].
 * \li "doInvert"  OBIT_bool (1,1,1) If TRUE invert solution [def FALSE];
 * \li "timeInt"   OBIT_float (1,1,1) Data integration time in sec [def 10 sec].
 * \param outUV    UV with which the output  Table is to be associated
 * \param ver      SN table version
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created ObitTableSN object which is 
 *                 associated with outUV.
 */
ObitTableSN* ObitTableSNGetZeroFR (ObitUV *inUV, ObitUV *outUV, olong ver, 
				   ObitErr *err)
{
  ObitTableSN *outCal=NULL;
  ObitTableSNRow *row=NULL;
  ObitAntennaList *AntList=NULL;
  ObitTableAN *ANTable=NULL;
  ObitSourceList *SouList=NULL;
  ObitTableSU *SUTable=NULL;
  ObitUVDesc *desc=NULL;
  ObitIOAccess access;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat *rec, solInt, t0, sumTime, cbase, lastTime=-1.0, lastSource=-1.0, lastFQID=-1.0;
  ofloat delTime, curSou=-1.0, curFQ=-1.0;
  olong i, ia, lrec, maxant;
  olong  nTime, SubA=0, ant1, ant2, lastSubA=-1;
  oint numPol, numIF, numOrb, numPCal;
  odouble DecR=0.0, RAR=0.0, ArrLong, cosdec=0.0;
  gboolean doCalSelect, doFirst=TRUE, someData=FALSE, gotAnt[MAXANT], invert=FALSE;
  ObitIOCode retCode;
  gchar *tname;
  gchar *routine = "ObitTableSNGetZeroFR";
 
   /* error checks */
  if (err->error) return outCal;
  g_assert (ObitUVIsA(inUV));
  desc = inUV->myDesc;

  /* Calibration/selection wanted? */ 
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

  /* Invert solutions? */
  ObitInfoListGetTest(inUV->info, "doInvert", &type, dim, &invert);

  /* open UV data to fully instantiate if not already open */
  if ((inUV->myStatus==OBIT_Inactive) || (inUV->myStatus==OBIT_Defined)) {
    retCode = ObitUVOpen (inUV, access, err);
    if (err->error) Obit_traceback_val (err, routine, inUV->name, outCal);
  }
  lrec = inUV->myDesc->lrec;
  t0 = -1.0e20;

  /* Create output */
  if (desc->jlocs>=0)  numPol = MIN (2, desc->inaxes[desc->jlocs]);
  else                 numPol = 1;
  if (desc->jlocif>=0) numIF = desc->inaxes[desc->jlocif];
  else                 numIF = 1;
  tname  = g_strconcat ("Calibration for: ",inUV->name, NULL);
  outCal = newObitTableSNValue(tname, (ObitData*)outUV, &ver, OBIT_IO_WriteOnly,  
			       numPol, numIF, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, outCal);

  /* Sorted? */
  if (desc->isort[0]=='T') {
    outCal->myDesc->sort[0] = outCal->TimeCol+1;
    outCal->myDesc->sort[1] = outCal->antNoCol+1;
  }

  /* Get parameters for calibration */
  /* "Solution interval" default 10 sec */
  solInt = 10.0;
  ObitInfoListGetTest(inUV->info, "solInt", &type, dim, (gpointer*)&solInt);
  solInt /= 86400.0;  /* to days */

  /* Source info - SU table or from header*/
  ver      = 1;
  access   = OBIT_IO_ReadOnly;
  SUTable = newObitTableSUValue ("SU table", (ObitData*)outUV, 
				 &ver, access, numIF, err);
  if (SUTable==NULL) {
    /* Get from header */
    ObitPrecessUVRaDecApp (desc, &RAR, &DecR);
    RAR      *= DG2RAD;
    DecR     *= DG2RAD;
    cosdec   = cos(DecR);
  } else { /* Use SU table */
    SouList = ObitTableSUGetList (SUTable, err);
    if (err->error) Obit_traceback_val (err, routine, outUV->name, outCal);
    
    /* Cleanup */
    SUTable = ObitTableSUUnref(SUTable);
  }

  /* Antenna info */
  ver      = 1;
  access   = OBIT_IO_ReadOnly;
  numOrb   = 0;
  numPCal  = 0;
  ANTable = newObitTableANValue ("AN table", (ObitData*)outUV, 
				 &ver, access, numOrb, numPCal, err);
  if (ANTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with AN table");
  AntList = ObitTableANGetList (ANTable, err);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, outCal);

  /* Cleanup */
  ANTable = ObitTableANUnref(ANTable);

  /* Antenna coordinates to wavelengths at reference frequency */
  for (i=0; i<AntList->number; i++) {
    AntList->ANlist[i]->AntXYZ[0] *= desc->crval[desc->jlocf]/VELIGHT;
    AntList->ANlist[i]->AntXYZ[1] *= desc->crval[desc->jlocf]/VELIGHT;
    AntList->ANlist[i]->AntXYZ[2] *= desc->crval[desc->jlocf]/VELIGHT;
  }
  ArrLong = AntList->ANlist[0]->AntLong;  /* Array longitude */

  /* Open table */
  if ((ObitTableSNOpen (outCal, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening input CL table", routine);
    return outCal;
  }

  /* Create Row */
  row = newObitTableSNRow (outCal);

  /* Attach row to output buffer */
  ObitTableSNSetRow (outCal, row, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outCal);

  /* Initialize */
  row->Time   = 0.0;
  row->TimeI  = 0.0;
  row->SourID = 0;
  row->antNo  = 0;
  row->SubA   = 0;
  row->NodeNo = 0;
  row->FreqID = 0;
  row->IFR       = 0.0;
  row->MBDelay1  = 0.0;
  /* IF dependent things */
  for (i=0; i<numIF; i++) {
    row->Real1[i]   = 1.0;
    row->Imag1[i]   = 0.0;
    row->Rate1[i]   = 0.0;
    row->Delay1[i]  = 0.0;
    row->Weight1[i] = 1.0;
    row->RefAnt1[i] = 0;
  }
  /* Multiple polarizations */
  if (numPol>1) {
    /* IF dependent things */
    for (i=0; i<numIF; i++) {
      row->Real2[i]   = 1.0;
      row->Imag2[i]   = 0.0;
      row->Rate2[i]   = 0.0;
      row->Delay2[i]  = 0.0;
      row->Weight2[i] = 1.0;
      row->RefAnt2[i] = 0;
    }
  } /* end two poln */
  /* List of antennas found */
  for (i=0; i<MAXANT; i++) gotAnt[i] = FALSE;
  maxant = AntList->number;

  /* Number of antennas */
  outCal->numAnt = maxant;

  /* Averaging time (sec) for time smearing correction */
  /* "Solution interval" default 10 sec */
  delTime = 10.0;
  ObitInfoListGetTest(inUV->info, "timeInt", &type, dim, &delTime);
  delTime /= 86400.0;  /* to days */

  /* loop looking at data */
  retCode = OBIT_IO_OK;
  sumTime = 0.0;
  nTime   = 0;
  doFirst = TRUE;
  while (retCode == OBIT_IO_OK) {
    
    /* read buffer */
    if (doCalSelect) retCode = ObitUVReadSelect (inUV, NULL, err);
    else retCode = ObitUVRead (inUV, NULL, err);
    if (err->error) Obit_traceback_val (err, routine, inUV->name, outCal);
    if (retCode==OBIT_IO_EOF) break; /* done? */
    
    /* Record pointer */
    rec = inUV->buffer;
    
    /* First time */
    if (t0<-1.0e10) {
      t0         = rec[inUV->myDesc->iloct];
      if (inUV->myDesc->ilocsu>0) lastSource = rec[inUV->myDesc->ilocsu];
      if (inUV->myDesc->ilocfq>0) lastFQID   = rec[inUV->myDesc->ilocfq];
      lastTime   = rec[inUV->myDesc->iloct];
      cbase      = rec[inUV->myDesc->ilocb]; /* Baseline */
      ant1       = (cbase / 256.0) + 0.001;
      ant2       = (cbase - ant1 * 256) + 0.001;
      lastSubA   = (olong)(100.0 * (cbase -  ant1 * 256 - ant2) + 0.5);
    }
    
    /* Loop over buffer */
    for (i=0; i<inUV->myDesc->numVisBuff; i++) {

      /* Accumulation or scan finished? If so, write "calibration".*/
      if (inUV->myDesc->ilocsu>0) curSou = rec[inUV->myDesc->ilocsu];
      if (inUV->myDesc->ilocfq>0) curFQ  = rec[inUV->myDesc->ilocfq];
      if ((rec[inUV->myDesc->iloct] > (t0+solInt)) || 
	  (curSou != lastSource) ||  (curFQ != lastFQID)) {
	
	/* Not first time - assume first descriptive parameter never blanked */
	if (nTime>0) {
	  /* if new scan write end of last scan and this time */
	  if ((curSou != lastSource) || (curFQ != lastFQID)) {
	    /* Need first entry for scan? */
	    if (doFirst) {
	      doFirst = FALSE;
	      row->Time  = t0;
	      row->TimeI = 0.0;
	      row->SourID = (oint)(MAX(lastSource,0.0)+0.5);
	      row->FreqID = (oint)(MAX(lastFQID,0.0)+0.5);
	      row->SubA   = lastSubA;
	      /* calculate/write rows */
	      SetRow (AntList, SouList, ArrLong, desc, &RAR, &DecR, delTime,
		      maxant,  gotAnt, invert, row, outCal, err);
	      if (err->error) Obit_traceback_val (err, routine, inUV->name, outCal);
	    } else { /* Not first scan */
	      /* values for end of previous scan */
	      row->Time   = lastTime; 
	      row->TimeI  = 0.0;
	      row->SourID = (oint)(MAX(lastSource,0.0)+0.5);
	      row->FreqID = (oint)(MAX(lastFQID,0.0)+0.5);
	      row->SubA   = lastSubA;
	      /* calculate/write rows */
	      SetRow (AntList, SouList, ArrLong, desc, &RAR, &DecR, delTime,
		      maxant,  gotAnt, invert, row, outCal, err);
	      if (err->error) Obit_traceback_val (err, routine, inUV->name, outCal);

	      /* Values for start of next scan */
	      row->Time   = rec[inUV->myDesc->iloct]; 
	      row->TimeI  = 0.0;
	      row->SourID = (oint)(rec[inUV->myDesc->ilocsu]+0.5);
	      row->SubA   = SubA;
	    } /* end write beginning of scan value */
	  } else {  /* in middle of scan - use average time */
	    /* Set descriptive info on Row */
	    row->Time  = sumTime/nTime;  /* time */
	    row->TimeI = 2.0 * (row->Time - t0);
	    row->SourID = (oint)(rec[inUV->myDesc->ilocsu]+0.5);
	    row->SubA   = SubA;
	  }
      
	  /* Write Cal table */
	  /* calculate/write rows */
	  SetRow (AntList, SouList, ArrLong, desc, &RAR, &DecR, delTime,
		  maxant,  gotAnt, invert, row, outCal, err);
	  if (err->error) Obit_traceback_val (err, routine, inUV->name, outCal);

	  /* initialize accumulators */
	  t0         = rec[inUV->myDesc->iloct];
	  if (inUV->myDesc->ilocsu>0) lastSource = rec[inUV->myDesc->ilocsu];
	  if (inUV->myDesc->ilocfq>0) lastFQID   = rec[inUV->myDesc->ilocfq];
	  sumTime    = 0.0;
	  nTime      = 0;
	  lastSubA   = -1;
 
	  /* Clear list of antennas found */
	  for (ia=0; ia<maxant; ia++) gotAnt[ia] = FALSE;

	} /* end of write entry if there is data */
      } /* end write entry */
      
      /* accumulate statistics
	 Antennas etc. */
      cbase = rec[inUV->myDesc->ilocb]; /* Baseline */
      ant1 = (cbase / 256.0) + 0.001;
      ant2 = (cbase - ant1 * 256) + 0.001;
      SubA = (olong)(100.0 * (cbase -  ant1 * 256 - ant2) + 0.5);
      if(lastSubA<=0) lastSubA = SubA;
      gotAnt[ant1] = TRUE;
      gotAnt[ant2] = TRUE;
      sumTime += rec[inUV->myDesc->iloct];
      lastTime = rec[inUV->myDesc->iloct];
      nTime++; /* how many data points */
      rec += inUV->myDesc->lrec; /* Data record pointer */
      someData = TRUE;
      
    } /* end loop over buffer load */
  } /* end loop reading data */
  
  /* Finish up any data in accumulator */
  if (nTime>0) {
    /* Set descriptive info on Row */
    row->Time   = sumTime/nTime;
    row->TimeI  = lastTime - t0;

    /* Write Cal table */
    /* calculate/write rows */
    SetRow (AntList, SouList, ArrLong, desc, &RAR, &DecR, delTime,
	    maxant,  gotAnt, invert, row, outCal, err);
    if (err->error) Obit_traceback_val (err, routine, inUV->name, outCal);
  } /* End final cal */

  /* Close cal table */
  if ((ObitTableSNClose (outCal, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing CL Table", routine);
    return outCal;
  }
  
  /* Close data */
  retCode = ObitUVClose (inUV, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outCal);

  /* Give warning if no data selected */
  if (!someData) Obit_log_error(err, OBIT_InfoWarn, 
				"%s: Warning: NO data selected", routine);
  /* Cleanup */
  AntList = ObitAntennaListUnref(AntList);
  SouList = ObitSourceListUnref(SouList);
  row     = ObitTableSNRowUnref(row); 
  return outCal;
} /* end ObitTableSNGetZeroFR */

/*---------------Private functions ----------------*/
/**
 * Sinc (sin pi*arg /pi*arg)
 * \param arg argument of function
 * \return sinc(arg).
 */
static ofloat Sinc (ofloat arg)
{
  ofloat out = 1.0;
  if (fabs(arg)<1.0e-5) return out;
  out = sin (G_PI*arg) / (G_PI*arg);
  return out;
} /* End Sinc */

/**
 * Calculate/write SN table entry for Zero Fringe Rate.
 * \param AntList Antenna list
 * \param SouList Source list or NULL
 * \param ArrLong Array longitude (rad)
 * \param desc    UV data descriptor
 * \param RAR     Apparent RA (rad)
 * \param DecR    Apparent declination (rad)
 * \param delTime UV data integrationt time (days)
 * \param maxant  Maximum antenna number
 * \param gotAnt  Array of flags for antennas with data
 * \param invert  If true invert solution
 * \param row     SN table row to use
 * \param outCal  Output SN table
 * \param err     ObitError/message stack
 */
static void SetRow (ObitAntennaList *AntList, ObitSourceList *SouList, 
		    odouble ArrLong, ObitUVDesc *desc, 
		    odouble *RAR, odouble *DecR, ofloat delTime,
		    olong maxant, gboolean  gotAnt[MAXANT], gboolean invert,
		    ObitTableSNRow *row, ObitTableSN *outCal, ObitErr *err) 
{
  odouble AntLst, HrAng=0.0, cosdec=0.0;
  olong ia, j, suid, iRow, numIF, numPol;
  gdouble twopi = 2.0* G_PI, omegaE=7.29115e-5;
  ofloat wt, uvw[3], bl[3], delay, phase;
  gchar *routine = "ObitTableSN:SetRow";

  numIF  = outCal->numIF; 
  numPol = outCal->numPol; 

  /* LST and hour angle (radians) */
  if (SouList) {  /* Source table? */
    suid     = row->SourID-1;
    *RAR     = SouList->SUlist[suid]->RAApp*DG2RAD;
    *DecR    = SouList->SUlist[suid]->DecApp*DG2RAD;
  }
  cosdec   = cos(*DecR);

  AntLst = AntList->GSTIAT0 + ArrLong + row->Time*AntList->RotRate;
  HrAng  = AntLst - *RAR;
  /* Loop over antennas found */
  for (ia=1; ia<=maxant; ia++) {
    if (!gotAnt[ia]) continue;
    bl[0] = (ofloat)AntList->ANlist[ia-1]->AntXYZ[0];
    bl[1] = (ofloat)AntList->ANlist[ia-1]->AntXYZ[1];
    bl[2] = (ofloat)AntList->ANlist[ia-1]->AntXYZ[2];
    /* Compute uvw - short baseline approximation */
    ObitUVUtilUVW (bl, *DecR, (ofloat)HrAng, uvw);
    /* Loop over IFs */
    for (j=0; j<numIF; j++) {
      /* IF freq */
      row->Rate1[j]   = +uvw[0]*cosdec*omegaE/desc->freq;
      row->Delay1[j]  = +uvw[2]/desc->freq;
      delay = uvw[2] * desc->freqIF[j]/desc->freq;
      phase = -twopi*(delay-(olong)delay);
      if (invert) { /* invert solution? */
	row->Rate1[j] = -row->Rate1[j];
	row->Delay1[j] = -row->Delay1[j];
	phase = -phase;
      }
      wt = fabs(Sinc(delTime*row->Rate1[j]) * Sinc(desc->chIncIF[j]*row->Delay1[j]));
      row->Real1[j]   = cos(phase);
      row->Imag1[j]   = sin(phase);
      row->Weight1[j] = wt;
      /* Multiple polarizations */
      if (numPol>1) {
	row->Real2[j]   = row->Real1[j];
	row->Imag2[j]   = row->Imag1[j];
	row->Rate2[j]   = row->Rate1[j];
	row->Delay2[j]  = row->Delay1[j];
	row->Weight2[j] = row->Weight1[j];
      }
      iRow = -1;
      row->antNo = ia;
      if ((ObitTableSNWriteRow (outCal, iRow, row, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "%s: ERROR writing SN Table file", routine);
	return ;
      }
    }
  } /* end loop over antennas */
} /* end SetRow */
