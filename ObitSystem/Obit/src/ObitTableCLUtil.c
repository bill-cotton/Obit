/* $Id$  */
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
#include "ObitTableCLUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableCLUtil.c
 * ObitTableCL class utility function definitions.
 */

/*----------------------Public functions---------------------------*/


/**
 * Copies CL tables from inUV to outUV with selection in inUV
 * If calibration is selected on inUV, no tables are copied
 * \param inUV     Input UV to copy from
 * \param outUV    Output UV to copy to
 * \param *err     ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableCLSelect (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableCL    *inTab=NULL, *outTab=NULL;
  ObitTableCLRow *inRow=NULL, *outRow=NULL;
  ObitInfoType type;
  olong itemp, iif, oif, i;
  olong highCLver, iCLver, inCLRow, outCLRow;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  union ObitInfoListEquiv InfoReal; 
  oint numPol, numIF, numTerm;
  gboolean wanted;
  gchar *CLType = "AIPS CL";
  gchar *routine = "ObitTableCLSelect";

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

  /* How many CL tables  */
  highCLver = ObitTableListGetHigh (inUV->tableList, CLType);

  /* Are there any? */
  if (highCLver <= 0) return OBIT_IO_OK;

  /* Loop over CL tables */
  for (iCLver=1; iCLver<=highCLver; iCLver++) {

    /* Get input table */
    numPol = 0;
    numIF  = 0;
    numTerm = 0;
    inTab = 
      newObitTableCLValue (inUV->name, (ObitData*)inUV, &iCLver, OBIT_IO_ReadOnly, 
			   numPol, numIF, numTerm, err);
    if (err->error) Obit_traceback_val (err, routine, inTab->name, retCode);
    /* Find it */
    if (inTab==NULL) continue;  /* No keep looping */

    /* Open input table */
    retCode = ObitTableCLOpen (inTab, OBIT_IO_ReadOnly, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inTab->name, retCode);

    /* Delete any old output table */
    retCode = ObitDataZapTable ((ObitData*)outUV, CLType, iCLver, err);
    if (err->error) Obit_traceback_val (err, routine, outUV->name, retCode);
 
    /* Create output table */
    numPol = MIN (2, inUV->mySel->numberPoln);
    numIF  = inUV->mySel->numberIF;
    numTerm = inTab->numTerm;
    outTab = 
      newObitTableCLValue (outUV->name, (ObitData*)outUV, &iCLver, OBIT_IO_WriteOnly, 
			   numPol, numIF, numTerm, err);
    if (err->error) Obit_traceback_val (err, routine, outUV->name, retCode);
    /* Create it? */
    Obit_retval_if_fail((outTab!=NULL), err, retCode,
			"%s: Could not create CL table %d for %s", 
			routine, iCLver, outTab->name);
 
  /* Open output table */
  retCode = ObitTableCLOpen (outTab, OBIT_IO_WriteOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inTab->name, retCode);
  
   /* Update header info */
    outTab->revision = inTab->revision;
    outTab->numAnt   = inTab->numAnt;
    outTab->mGMod    = inTab->mGMod;

    /* Set rows */
    inRow  = newObitTableCLRow (inTab);
    outRow = newObitTableCLRow (outTab);
    ObitTableCLSetRow (outTab, outRow, err);
    if (err->error) Obit_traceback_val (err, routine, outTab->name, retCode);

    /* Loop over table copying selected data */
    outCLRow = -1;
    for (inCLRow=1; inCLRow<=inTab->myDesc->nrow; inCLRow++) {
      retCode = ObitTableCLReadRow (inTab, inCLRow, inRow, err);
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
      outRow->atmos   = inRow->atmos;
      outRow->Datmos   = inRow->Datmos;
      for (i=0; i<inTab->numTerm; i++) 
	outRow->GeoDelay[i] = inRow->GeoDelay[i];
      outRow->MBDelay1  = inRow->MBDelay1;
      outRow->clock1    = inRow->clock1;
      outRow->Dclock1   = inRow->Dclock1;
      outRow->dispers1  = inRow->dispers1;
      outRow->Ddispers1 = inRow->Ddispers1;
      oif = 0;
      for (iif=inUV->mySel->startIF-1; 
	   iif<inUV->mySel->startIF+inUV->mySel->numberIF-1;
	   iif++) {
	     outRow-> DopplerOff[oif] = inRow-> DopplerOff[iif];
	     outRow->Real1[oif]   = inRow->Real1[iif];
	     outRow->Imag1[oif]   = inRow->Imag1[iif];
	     outRow->Delay1[oif]  = inRow->Delay1[iif];
	     outRow->Rate1[oif]   = inRow->Rate1[iif];
	     outRow->Weight1[oif] = inRow->Weight1[iif];
	     outRow->RefAnt1[oif] = inRow->RefAnt1[iif];
	     oif++;
	   }
      if (numPol>1) {
	outRow->MBDelay2  = inRow->MBDelay2;
	outRow->clock2    = inRow->clock2;
	outRow->Dclock2   = inRow->Dclock2;
	outRow->dispers2  = inRow->dispers2;
	outRow->Ddispers2 = inRow->Ddispers2;
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

      retCode = ObitTableCLWriteRow (outTab, outCLRow, outRow, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, inUV->name, retCode);
    } /* end loop over rows */
    
    /* Close tables */
    retCode = ObitTableCLClose (inTab, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inTab->name, retCode);
    retCode = ObitTableCLClose (outTab, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, outTab->name, retCode);
 
    /* release table objects */
    inTab  = ObitTableCLUnref(inTab);
    outTab = ObitTableCLUnref(outTab);

    /* release row objects */
    inRow  = ObitTableCLRowUnref(inRow);
    outRow = ObitTableCLRowUnref(outRow);
  } /* end loop over tables */

  return retCode;
} /* end ObitTableCLSelect */

