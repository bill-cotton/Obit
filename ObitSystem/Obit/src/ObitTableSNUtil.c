/* $Id: ObitTableSNUtil.c,v 1.5 2007/07/07 16:21:09 bcotton Exp $  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005                                               */
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

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableSNUtil.c
 * ObitTableSN class utility function definitions.
 */

/*----------------------Public functions---------------------------*/


/**
 * Copies SN tables from inUV to outUV with selection in inUV
 * If calibration is selected on inUV, no tables are copied
 * \param inUV     Input UV to copy from
 * \param outUV    Output UV to copy to
 * \param *err     ObitErr error stack.
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
 * \return pointer to the new object.
 */
ObitTableSN* ObitTableSNUtilInvert (ObitTableSN *inSN, ObitData *outData, olong *outVer, 
				    ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  ObitTableSN    *outSN=NULL;
  ObitTableSNRow *inRow=NULL, *outRow=NULL;
  ofloat amp, phase, fblank = ObitMagicF();
  olong inSNRow, outSNRow, oif, iif;
  gchar *routine = "ObitTableSNUtilInvert";
  
  /* error checks */
  if (err->error) return outSN;
  g_assert (ObitTableSNIsA(inSN));
  g_assert (ObitDataIsA(outData));
  
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
	if ((inRow->Real1[iif]!=fblank) && (inRow->Imag1[iif]!=fblank)) {
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
	} else { /* blanked */
	  outRow->Real1[oif]   = fblank;
	  outRow->Imag1[oif]   = fblank;
	  outRow->Delay1[oif]  = fblank;
	  outRow->Rate1[oif]   = fblank;
	}
	outRow->Weight1[oif] = inRow->Weight1[iif];
	outRow->RefAnt1[oif] = inRow->RefAnt1[iif];
	oif++;
      }
      if (inSN->numPol>1) {
	outRow->MBDelay2 = -inRow->MBDelay2;
	oif = 0;
	for (iif=0; iif<inSN->numIF; iif++) {
	  if ((inRow->Real2[iif]!=fblank) && (inRow->Imag2[iif]!=fblank)) {
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
	  } else { /* blanked */
	    outRow->Real2[oif]   = fblank;
	    outRow->Imag2[oif]   = fblank;
	    outRow->Delay2[oif]  = fblank;
	    outRow->Rate2[oif]   = fblank;
	  }
	  outRow->Weight2[oif] = inRow->Weight2[iif];
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
