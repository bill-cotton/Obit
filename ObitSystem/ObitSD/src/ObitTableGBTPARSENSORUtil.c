/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009                                               */
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

#include "ObitTableGBTPARSENSORUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableGBTPARSENSORUtil.c
 * ObitTableGBTPARSENSOR class utility function definitions.
 */

/*----------------------Public functions---------------------------*/


/**
 * Determine Get translation table between index in dac data vectors and the
 * table of detector position offsets.
 * \param in     GBTPARSENSOR table
 * \param n      Length of row, col
 * \param row    row numbers (0-rel) in order of detector offset table
 * \param col    column numbers (0-rel) in order of detector offset table
 * \param *err   ObitErr error stack.
 * \return array of olong of length n, entries are the (0-rel) indices in the 
 *         detector offset table corresponding to each value in the dac arrays
 *         -1 -> not found.
 */
olong* ObitTableGBTPARSENSORXlate (ObitTableGBTPARSENSOR *in, olong n,
				   olong *row, olong *col, ObitErr *err)
{
  olong *out=NULL;
  ObitIOCode retCode;
  ObitTableGBTPARSENSORRow *Row = NULL;
  olong i=0, j, irow, *tcols=NULL, *trows=NULL;
  gchar *routine = "ObitTableGBTPARSENSORXlate";

  /* error checks */
  if (err->error) return out;

  /* Open table */
  retCode = ObitTableGBTPARSENSOROpen (in, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine,in->name, out);
 
  /* allocate arrays */
  out   = g_malloc0(MAX(n,in->myDesc->nrow)*sizeof(olong));
  trows = g_malloc0(in->myDesc->nrow*sizeof(olong));
  tcols = g_malloc0(in->myDesc->nrow*sizeof(olong));

  /* Create table row */
  Row  = newObitTableGBTPARSENSORRow (in);

  /* Swallow table */
  /* Loop over table checking times */
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTableGBTPARSENSORReadRow (in, irow, Row, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
    trows[irow-1] = Row->row;
    tcols[irow-1] = Row->col;    
  } /* end loop over rows */

  /* Close table */
  retCode = ObitTableGBTPARSENSORClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) {
    Obit_traceback_val (err, routine, in->name, out);
  }

  /* release row object */
  Row = ObitTableGBTPARSENSORRowUnref(Row);

  /* Find lookup table  */
  for (i=0; i<MAX(n,in->myDesc->nrow); i++) out[i] = -1;  /* init */
  for (i=0; i<n; i++) {                   /* Loop over input */
    for (j=0; j<in->myDesc->nrow; j++) {  /* Loop over table */
      /* match? */
      if ((row[i]==trows[j]) && (col[i]==tcols[j])) {
	out[j] = i;
	break;
      }
    }
  }

  if (trows) g_free(trows);
  if (trows) g_free(tcols);

  return out;
} /* end ObitTableGBTPARSENSORWantSour */

