/* $Id: ObitTableCCUtil.c 509 2015-03-25 19:53:54Z bill.cotton $   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2016                                               */
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
#include "ObitTableBPUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableBPUtil.c
 * ObitTableBP class utility function definitions.
 */

/*----------------------Private function prototypes----------------------*/
/*----------------------Public functions---------------------------*/
/**
 * Copy the contents of table in to the end of out.
 * Drops deselected entries (both AIPS and FITS).
 * \param in     Input BP table
 * \param out    Output BP table
 * \param *err   ObitErr error stack.
 */
void ObitTableBPUtilAppend (ObitTableBP *in, ObitTableBP *out, ObitErr *err)
{
  ObitTableBPRow *row = NULL;
  olong irow, orow;
  gchar *routine = "ObitTableBPUtilAppend";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableBPIsA(in), err, 
		       "%s input %s not an BP Table", routine, in->name);
  Obit_return_if_fail (ObitTableBPIsA(out), err, 
		       "%s output %s not an BP Table", routine, out->name);
  
  /* Open output */
  ObitTableBPOpen (out, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);
  row  = newObitTableBPRow (out);
  ObitTableSetRow ((ObitTable*)out, (ObitTableRow*)row, err);

  /* Will not be sorted */
  out->myDesc->sort[0] = 0;
  out->myDesc->sort[1] = 0;

  /* Open input */
  ObitTableBPOpen (in, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Loop over input Table */
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTableBPReadRow (in, irow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (row->status<0) continue;  /* Skip deselected record */

    /* Write row */
    orow = -1;
    ObitTableBPWriteRow (out, orow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, out->name);
  } /* End loop over table */
    
  /* cleanup/close up */
  row = ObitTableBPRowUnref(row);
  ObitTableBPClose (in,  err);
  ObitTableBPClose (out, err);
} /* end ObitTableBPUtilAppend */
