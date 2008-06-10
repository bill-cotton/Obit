/* $Id: ObitTableNXUtil.c,v 1.3 2007/08/17 01:10:42 bcotton Exp $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006                                               */
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

#include "ObitTableNXUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableNXUtil.c
 * ObitTableNX class utility function definitions.
 */

/*----------------------Public functions---------------------------*/


/**
 * Determine if a given source ID has data in a given timerange.
 * \param in        NX table to access , opened and closed externally
 * \param sid       Desired Source id
 * \param timerange Timerange in question
 * \param *err      ObitErr error stack.
 * \return TRUE if iNdeX table shows data for the given source 
 *   in the given timerange, else FALSE
 */
gboolean ObitTableNXWantSour (ObitTableNX *in, olong sid, 
			      ofloat timerange[2], ObitErr *err)
{
  gboolean out = FALSE;
  ObitTableNXRow *row = NULL;
  gboolean want;
  olong irow;
  gchar *routine = "ObitTableNXWantSour ";

  /* error checks */
  if (err->error) return out;

  row  = newObitTableNXRow (in);

  /* Loop over table checking times */
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTableNXReadRow (in, irow, row, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
    if (row->status==-1) continue;  /* Deselected? */
    
    /* Want this one? */
    want = (sid==row->SourID) &&
      /* Is either beginning or end of timerange in scan */
      (((timerange[0]>=(row->Time-0.5*row->TimeI)) &&
	((timerange[0]<=(row->Time+0.5*row->TimeI)))) ||
       ((timerange[1]>=(row->Time-0.5*row->TimeI)) &&
	((timerange[1]<=(row->Time+0.5*row->TimeI)))) ||
       /* Or scan entirely inside timerange */
       ((timerange[0]<=(row->Time-0.5*row->TimeI)) &&
	(timerange[1]>=(row->Time+0.5*row->TimeI))))
       ;
    
    if (want) {out = TRUE; break;}
  } /* end loop over rows */

  /* release row object */
  row = ObitTableNXRowUnref(row);

  return out;
} /* end ObitTableNXWantSour */

