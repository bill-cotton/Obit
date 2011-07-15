/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2011                                          */
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

#include "ObitTableFQUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableFQUtil.c
 * ObitTableFQ class utility function definitions.
 */

/*----------------------Public functions---------------------------*/


/**
 * Reads and returns information for a specified FQid from "AIPS FQ" table..
 * \param in       FQ table to access 
 * \param fqid     Desired frequency id 
 * \param nif      Number of IFs in data  
 * \param freqOff  IF Frequency offsets from reference frequency (Hz)
 *                 returns pointer to newly malloced memory, must be g_freeed.
 * \param sideBan  Sideband array, 1=> upper, -1=lower (?)
 *                 returns pointer to newly malloced memory, must be g_freeed.
 * \param chBandw  channel bandwidths per IF (Hz),
 *                 returns pointer to newly malloced memory, must be g_freeed.
 * \param *err     ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableFQGetInfo (ObitTableFQ *in, oint fqid, oint *nif,
			       odouble **freqOff, oint **sideBand,
			       ofloat **chBandw, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong nRowPIO;
  odouble *freqOffOut, *dRow;
  oint *sideBandOut, i, *iRow;
  ofloat *chBandwOut, *fRow;
  gchar *routine = "ObitTableFQGetInfo";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitTableFQIsA(in));

  /* Read one row  */
  nRowPIO = 1;
  ObitInfoListPut(in->info, "nRowPIO", OBIT_long, dim, 
		  (gconstpointer)&nRowPIO, err);
  if (err->error) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Open */
  retCode = ObitTableFQOpen (in, OBIT_IO_ReadOnly,  err);
  if (err->error) /* add traceback,return */
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Number of IFs */
  *nif = in->numIF;
  
  /* read row fqid */
  retCode = ObitTableRead ((ObitTable*)in, fqid, NULL, err);
  if (err->error) 
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Copy data */
  dRow = (odouble*)in->buffer;
  iRow = (oint*)in->buffer;
  fRow = (ofloat*)in->buffer;
  freqOffOut  = g_malloc0(in->myDesc->repeat[in->freqOffCol]*sizeof(odouble));
  for (i=0; i<in->myDesc->repeat[in->freqOffCol]; i++) 
    freqOffOut[i] = dRow[in->freqOffOff+i];
  
  sideBandOut = g_malloc0(in->myDesc->repeat[in->sideBandCol]*sizeof(olong));
  for (i=0; i<in->myDesc->repeat[in->sideBandCol]; i++) 
    sideBandOut[i] = iRow[in->sideBandOff+i];
  
  chBandwOut  = g_malloc0(in->myDesc->repeat[in->chWidthCol]*sizeof(ofloat));
  for (i=0; i<in->myDesc->repeat[in->chWidthCol]; i++) 
    chBandwOut[i] = fRow[in->chWidthOff+i];
  
  /* Set output pointers */ 
  *freqOff  = freqOffOut;
  *sideBand = sideBandOut;
  *chBandw  = chBandwOut;
  
  /* Close */
  retCode = ObitTableFQClose (in,  err);
  if (err->error) 
    Obit_traceback_val (err, routine, in->name, retCode);
  
  return retCode;
} /* end ObitTableFQGetInfo */

/**
 * Writes information for a specified FQid to "AIPS FQ" table.
 * \param in       FQ table to access
 * \param fqid     Desired frequency id
 * \param nif      Number of IFs in data
 * \param freqOff  IF Frequency offsets from reference frequency (Hz)
 * \param sideBan  Sideband array, 1=> upper, -1=lower (?)
 * \param chBandw  channel bandwidths per IF (Hz),
 * \param *err     ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableFQPutInfo (ObitTableFQ *in, oint fqid, oint nif,
			       odouble *freqOff, oint *sideBand,
			       ofloat *chBandw, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
   gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
   olong nRowPIO, *siRow;
   odouble *dRow;
   oint   *iRow, oitemp;
   ofloat *fRow;
   gchar *routine = "ObitTableFQPutInfo";
   
   /* error checks */
   g_assert (ObitErrIsA(err));
   if (err->error) return retCode;
   g_assert (ObitTableFQIsA(in));

   /* Number of IFs */
   in->numIF = nif;

   /* Open */
   retCode = ObitTableFQOpen (in, OBIT_IO_ReadWrite,  err);
   if (err->error) 
     Obit_traceback_val (err, routine, in->name, retCode);
   
   /* Save on info */
   oitemp = (oint)in->numIF;
   ObitInfoListAlwaysPut(in->myDesc->info, "NO_IF", OBIT_oint, dim, 
			 (gconstpointer)&oitemp);
   
   /* Write one row  */
   nRowPIO = 1;
   ObitInfoListPut(in->info, "nRowPIO", OBIT_long, dim, 
		   (gconstpointer)&nRowPIO, err);
   if (err->error) /* add traceback,return */
     Obit_traceback_val (err, routine, in->name, retCode);
   
   /* read row fqid */
   retCode = ObitTableRead ((ObitTable*)in, fqid, NULL,  err);
   if (err->error) 
     Obit_traceback_val (err, routine, in->name, retCode);
   
   /* Copy data */
   dRow = (odouble*)in->buffer;
   siRow = (olong*)in->buffer;
   iRow = (oint*)in->buffer;
   fRow = (ofloat*)in->buffer;
   for (i=0; i<in->myDesc->repeat[in->freqOffCol]; i++) 
     dRow[in->freqOffOff+i] = freqOff[i];
   
   for (i=0; i<in->myDesc->repeat[in->sideBandCol]; i++) 
     iRow[in->sideBandOff+i] = sideBand[i];
   
   for (i=0; i<in->myDesc->repeat[in->chWidthCol]; i++) 
     fRow[in->chWidthOff+i] = chBandw[i];
   
   /* Mark as modified */
   siRow[in->myDesc->statusOff] = 1;
   in->myDesc->numRowBuff = 1;
   
   /* rewrite */
   retCode = ObitTableWrite ((ObitTable*)in, fqid, NULL,  err);
   if (err->error) 
     Obit_traceback_val (err, routine, in->name, retCode);
   
   /* Close */
   retCode = ObitTableFQClose (in,  err);
   if (err->error) 
     Obit_traceback_val (err, routine, in->name, retCode);
   
   return retCode;
} /* end ObitTableFQPutInfo */


/**
 * Copies FQ table from inUV to outUV with selection in inUV
 * \param inUV     Input UV to copy from
 * \param outUV    Output UV to copy to
 * \param SouIFOff if NonNULL, source dependent values to be added to the 
 *                 IF frequencies. Includes selection by IF.
 * \param *err     ObitErr error stack.
 * \param SouBW    if >0.0, source dependent total bandwidth
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableFQSelect (ObitUV *inUV, ObitUV *outUV, odouble *SouIFOff,
			      odouble SouBW, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableFQ    *inTab=NULL, *outTab=NULL;
  ObitTableFQRow *inRow=NULL, *outRow=NULL;
  ObitErr *terr=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong iif, oif, nif, nchAvg, chInc, IFInc;
  olong iFQver, inFQRow, outFQRow, highFQver, maxIF=1;
  oint numIF;
  gboolean wanted;
  gchar *FQType = "AIPS FQ";
  gchar *routine = "ObitTableFQSelect";

  /* error checks */
  if (err->error) return retCode;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));

  /* Fully instantiate UV files */
  ObitUVFullInstantiate (inUV, TRUE, err);
  if (err->error )Obit_traceback_val (err, routine, inUV->name, retCode);
  /* Ignore any problems */
  ObitErrLog(err); /* Show any pending messages as they may get lost */
  ObitUVFullInstantiate (outUV, FALSE, err);
  ObitErrClear(err); 

  /* How many FQ tables  */
  highFQver = ObitTableListGetHigh (inUV->tableList, FQType);

  /* Are there any? */
  if (highFQver <= 0) return OBIT_IO_OK;

  /* How many channels to average */
  nchAvg = 1;
  chInc  = MAX (1, inUV->mySel->channInc);
  ObitInfoListGetTest(inUV->info, "NumChAvg",  &type, dim, &nchAvg);  
  nchAvg /= chInc;   /* Channel increment */
  nchAvg = MAX (1,nchAvg);
  nchAvg = MIN (nchAvg, inUV->myDesc->inaxes[inUV->myDesc->jlocf]);

  /* Should only be one FQ table */
  iFQver = 1;
  if (inUV->myDesc->jlocif>=0) {
    nif   = inUV->myDesc->inaxes[inUV->myDesc->jlocif];
    IFInc = MAX (1, inUV->mySel->IFInc);   /* IF increment */
    maxIF = inUV->mySel->startIF+inUV->mySel->numberIF*IFInc-1;
    maxIF = MIN (maxIF, 
		 ((ObitUVDesc*)inUV->myIO->myDesc)->inaxes[((ObitUVDesc*)inUV->myIO->myDesc)->jlocif]);
  } else {
    nif = 1;
    IFInc = 1;
  }

  /* Get input table */
  numIF = 0;
  inTab = 
    newObitTableFQValue (inUV->name, (ObitData*)inUV, &iFQver, OBIT_IO_ReadOnly, 
			 numIF, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, retCode);
  /* Find it? */
   Obit_retval_if_fail(((inTab!=NULL) || (nif<=1)), err, retCode,
		      "%s: Could not find FQ table for %s %d IFs", 
		      routine, inUV->name, nif);
  
  /* Open input table */
  retCode = ObitTableFQOpen (inTab, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inTab->name, retCode);
  
  /* Delete any old output table - ignore any problems */
  terr = newObitErr();
  retCode = ObitDataZapTable ((ObitData*)outUV, FQType, iFQver, terr);
  terr = ObitErrUnref(terr);
  
  /* Create output table */
  numIF  = inUV->mySel->numberIF;
  outTab = 
    newObitTableFQValue (outUV->name, (ObitData*)outUV, &iFQver, OBIT_IO_WriteOnly, 
			 numIF, err);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, retCode);
  /* Create it? */
  Obit_retval_if_fail((outTab!=NULL), err, retCode,
		      "%s: Could not create FQ table %d for %s", 
		      routine, iFQver, outTab->name);
  
  /* Open output table */
  retCode = ObitTableFQOpen (outTab, OBIT_IO_WriteOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inTab->name, retCode);
  
  /* Set rows */
  inRow  = newObitTableFQRow (inTab);
  outRow = newObitTableFQRow (outTab);
  ObitTableFQSetRow (outTab, outRow, err);
  if (err->error) Obit_traceback_val (err, routine, outTab->name, retCode);
  
  /* Loop over table copying selected data */
  outFQRow = -1;
  for (inFQRow=1; inFQRow<=inTab->myDesc->nrow; inFQRow++) {
    retCode = ObitTableFQReadRow (inTab, inFQRow, inRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inUV->name, retCode);
    if (inRow->status==-1) continue;
    
    /* Want this one? Use all if none specified */
    wanted = ((inUV->mySel->FreqID<=0) || 
	      (inRow->fqid==inUV->mySel->FreqID));
    if (!wanted) continue;
    
    /* Copy selected data */
    outRow->fqid     = inRow->fqid;
    oif = 0;
    for (iif=inUV->mySel->startIF-1; iif<maxIF; iif+=IFInc) {
      outRow->freqOff[oif]  = inRow->freqOff[iif] - 
	inRow->freqOff[inUV->mySel->startIF-1]; /* New reference freq */
      outRow->chWidth[oif]  = inRow->chWidth[iif] * nchAvg;
      outRow->totBW[oif]    = inRow->totBW[iif];
      outRow->sideBand[oif] = inRow->sideBand[iif];
      /* Source dependent OFFSets */
      if (SouIFOff!=NULL) outRow->freqOff[oif] += SouIFOff[oif];
      if (SouBW>0.0) outRow->totBW[oif] = SouBW;
      oif++;
    }
    
    retCode = ObitTableFQWriteRow (outTab, outFQRow, outRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inUV->name, retCode);
  } /* end loop over rows */
  
    /* Close tables */
  retCode = ObitTableFQClose (inTab, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inTab->name, retCode);
  retCode = ObitTableFQClose (outTab, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, outTab->name, retCode);
  
  /* release table objects */
  inTab  = ObitTableFQUnref(inTab);
  outTab = ObitTableFQUnref(outTab);
  
  /* release row objects */
  inRow  = ObitTableFQRowUnref(inRow);
  outRow = ObitTableFQRowUnref(outRow);

  return retCode;
} /* end ObitTableFQSelect */

