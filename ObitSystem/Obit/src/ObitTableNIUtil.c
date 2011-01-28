/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2008                                               */
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

#include "ObitTableNIUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableNIUtil.c
 * ObitTableNI class utility function definitions.
 */

/*----------------------Private prototypes---------------------------*/
/*----------------------Public functions---------------------------*/
/**
 * Swallow an NI table to an ObitFArray 
 * \param in Table to swallow
 * \param err Obit error stack object.
 * \return FArray with elaments (all floats):
 * \li [0] Time    - The center time (days)
 * \li [1] TimeI   - Time interval of the solution 
 * \li [2] antNo   - Antenna number, 0=$>$ al
 * \li [3] SourId  - Source number, 0=$>$ all 
 * \li [4] SubA    - Subarray number, 0=$>$ all 
 * \li [5] weight  - Weight
 * \li [6...] coef - Zernike model coeffients
 *    The number of coefficients is lrow-6.
 */
ObitFArray* ObitTableNIUtilSwallow (ObitTableNI *in, ObitErr *err)
{
  ObitFArray *out = NULL;
  ObitIOCode retCode;
  ObitTableNIRow *NIRow = NULL;
  olong nrow, irow, lrow, ncoef, icoef, naxis[2], pos[2];
  ofloat *Row;
  gchar *routine = "ObitTableNIUtilSwallow";
 
  /* error checks */
  if (err->error) return out;

  /* Open */
  retCode = ObitTableNIOpen (in, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);
  
  /* Create output structure */
  nrow  = in->myDesc->nrow;
  ncoef = in->numCoef;
  lrow  = 6 + ncoef;
  naxis[0] = lrow; naxis[1] = nrow;
  out = ObitFArrayCreate ("NI Table", 2, naxis);

  /* Create table row */
  NIRow = newObitTableNIRow (in);

  /* Loop over table reading CCs */
  pos[0] = 0;
  for (irow=1; irow<=nrow; irow++) {

    retCode = ObitTableNIReadRow (in, irow, NIRow, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
    if (NIRow->status<0) continue;   /* Valid row? */

    /* Save data */
    pos[1] = irow-1;
    Row = ObitFArrayIndex (out, pos);
    Row[0] = (ofloat)NIRow->Time;
    Row[1] =         NIRow->TimeI;
    Row[2] = (ofloat)NIRow->antNo;
    Row[3] = (ofloat)NIRow->SourId;
    Row[4] = (ofloat)NIRow->SubA;
    Row[5] =         NIRow->weight;
    for (icoef=0; icoef<ncoef; icoef++) Row[6+icoef] = NIRow->coef[icoef];
  } /* end loop over table */

  /* Close */
  retCode = ObitTableNIClose (in, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  /* Cleanup */
  NIRow = ObitTableNIRowUnref(NIRow);
  return out;
} /* end ObitTableNIUtilSwallow */

/**
 * Return  0-rel index of row in NI Table FArray prior to a given time
 * \param NIArray FArray with NI table returned by ObitTableNIUtilSwallow
 * Rows each contain (all floats):
 * \li [0] Time    - The center time (days)
 * \li [1] TimeI   - Time interval of the solution 
 * \li [2] antNo   - Antenna number, 0=$>$ al
 * \li [3] SourId  - Source number, 0=$>$ all 
 * \li [4] SubA    - Subarray number, 0=$>$ all 
 * \li [5] weight  - Weight
 * \li [6...] coef - Zernike model coeffients
 * \param time Target time in days, if time is prior to the first
 *    time in NIArray then returns index of the first row.
 *  \param prior start search at index prior
 * \return row index in NIArray
 */
olong ObitTableNIUtilPrior (ObitFArray* NIArray, ofloat time, olong *prior)
{
  olong out=*prior;
  olong nrow, irow, lrow, l1, pos[2]={0,0};
  ofloat *Row;

  /* How long are rows? */
  lrow = NIArray->naxis[0];

  /* How many rows? */
  nrow = NIArray->naxis[1];

  /* Row pointer */
  l1 = MAX (0, *prior);
  pos[1] = l1;
  Row = ObitFArrayIndex (NIArray, pos);
  out = l1;
  for (irow=l1; irow<nrow; irow++) {
    if (Row[0]>=time) break; /* Want one just before later time */
    out = irow;
    Row += lrow;
  }

  return out;
} /* end ObitTableNIUtilPrior */

/**
 * Return 0-rel index of row in NI Table FArray following a given time 
 * \param NIArray FArray with NI table returned by ObitTableNIUtilSwallow
 * Rows each contain (all floats):
 * \li [0] Time    - The center time (days)
 * \li [1] TimeI   - Time interval of the solution 
 * \li [2] antNo   - Antenna number, 0=$>$ al
 * \li [3] SourId  - Source number, 0=$>$ all 
 * \li [4] SubA    - Subarray number, 0=$>$ all 
 * \li [5] weight  - Weight
 * \li [6...] coef - Zernike model coeffients
 * \param time Target time in days, if time follows the last
 *    time in NIArray then returns index of the last row.
 *  \param follow start search at index follow
 * \return row index in NIArray
 */
olong ObitTableNIUtilFollow (ObitFArray* NIArray, ofloat time, olong *follow)
{
  olong out=*follow;
  olong nrow, irow, lrow, l1, pos[2]={0,0};
  ofloat *Row;

  /* How long are rows? */
  lrow = NIArray->naxis[0];

  /* How many rows? */
  nrow = NIArray->naxis[1];

  /* Row pointer */
  l1 = MAX (0, *follow);
  pos[1] = l1;
  Row = ObitFArrayIndex (NIArray, pos);
  out = l1;
  for (irow=l1; irow<nrow; irow++) {
    out = irow;
    if (Row[0]>time) break;
    Row += lrow;
  }

  return out;
} /* end ObitTableNIUtilFollow */

