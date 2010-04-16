/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006,2010                                          */
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

/**
 * Make flagging entries in FGTab for beginning and/or end of scans
 * for selected data.
 * \param in        NX table to use
 * Control parameter on info element of in:
 * \li "begDrop"   OBIT_float  (1,1,1) Time (min) at start of scan to flag
 * \li "endDrop"   OBIT_float  (1,1,1) Time (min) at end of scan to flag
 * \li "Reason"    OBIT_string (24,1,1) Reason String for FG table
 * \param FGTab     FG Table to write
 * \param sel       UV data selector to determine which data wanted.
 * \param maxAnt    Maximum antenna number
 * \param err       ObitErr error stack.
 */
void ObitTableNXUtilQuack (ObitTableNX *in, ObitTableFG *FGTab,
			   ObitUVSel *sel, olong maxAnt, ObitErr *err)
{
  ObitTableNXRow *NXrow = NULL;
  ObitTableFGRow *FGrow = NULL;
  gboolean want, *AntOK;
  olong irow, orow, iant, nscan=0;
  ofloat begDrop, endDrop, sumTime=0.0;
  gchar Reason[64];
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gchar *routine = "ObitTableNXUtilQuack ";

  /* error checks */
  if (err->error) return;

  /* Parameters */
  begDrop = 0.0;
  ObitInfoListGetTest(in->info, "begDrop", &type, dim, &begDrop);
  endDrop = 0.0;
  ObitInfoListGetTest(in->info, "endDrop", &type, dim, &endDrop);
  sprintf (Reason, "                        ");
  ObitInfoListGetTest(in->info, "Reason",  &type, dim, Reason);

  /* Anything to do? */
  if ((begDrop<=0.0) && (endDrop<=0.0)) return;

  /* to days */
  begDrop /= 1440.0;
  endDrop /= 1440.0;

  /* Which antennas wanted */
  AntOK  = g_malloc0(maxAnt*sizeof(gboolean));
  for (iant = 0; iant<maxAnt; iant++)
    AntOK[iant] = ObitUVSelWantAnt (sel, iant+1);
 
  /* Open Index table */
  if ((ObitTableNXOpen (in, OBIT_IO_ReadOnly, err)
       != OBIT_IO_OK) || (err->error)) goto cleanup;
  NXrow  = newObitTableNXRow (in);

  /* Open Flag table */
  if ((ObitTableFGOpen (FGTab, OBIT_IO_ReadWrite, err)
       != OBIT_IO_OK) || (err->error)) goto cleanup;
  FGrow  = newObitTableFGRow (FGTab);

  /* Attach  row to output buffer */
  ObitTableFGSetRow (FGTab, FGrow, err);
  if (err->error) goto cleanup;

  /* If there are entries in the table, mark it unsorted */
  if (FGTab->myDesc->nrow>0) 
    {FGTab->myDesc->sort[0]=0; FGTab->myDesc->sort[1]=0;}
  
  /* Init  Flag row */
  FGrow->SourID       = 0;
  FGrow->SubA         = 0;
  FGrow->freqID       = 0;
  FGrow->TimeRange[0] = 0;
  FGrow->TimeRange[1] = 0;
  FGrow->ants[0]      = 1;
  FGrow->ants[1]      = 0;
  FGrow->chans[0]     = sel->startChann;
  FGrow->chans[1]     = sel->startChann + sel->numberChann - 1;
  FGrow->ifs[0]       = sel->startIF;
  FGrow->ifs[1]       = sel->startIF + sel->numberIF - 1;
  FGrow->pFlags[0]    = 0;
  strncpy (FGrow->reason, Reason, 24);
  if ((sel->Stokes[0]==' ') || (sel->Stokes[0]==0)){
    FGrow->pFlags[0] = 15;
  } else {
    if ((sel->Stokes[0]=='R')&&(sel->Stokes[1]=='R')) FGrow->pFlags[0] += 1;
    if ((sel->Stokes[1]=='L')&&(sel->Stokes[1]=='L')) FGrow->pFlags[0] += 2;
    if ((sel->Stokes[2]=='R')&&(sel->Stokes[1]=='L')) FGrow->pFlags[0] += 4;
    if ((sel->Stokes[4]=='L')&&(sel->Stokes[1]=='R')) FGrow->pFlags[0] += 8;
    if (sel->Stokes[0]=='1') FGrow->pFlags[0] += 1;
    if (sel->Stokes[1]=='1') FGrow->pFlags[0] += 2;
    if (sel->Stokes[2]=='1') FGrow->pFlags[0] += 4;
    if (sel->Stokes[3]=='1') FGrow->pFlags[0] += 8;
  }
  
  /* Loop over NX table  */
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTableNXReadRow (in, irow, NXrow, err);
    if (err->error) goto cleanup;
    if (NXrow->status==-1) continue;  /* Deselected? */
    
    /* Want this one? */
    want = ObitUVSelWantSour (sel, NXrow->SourID);
    want = want && ((sel->SubA <= 0) || (sel->SubA ==  NXrow->SubA));
    want = want && ((sel->FreqID <= 0) || (sel->FreqID ==  NXrow->FreqID));
    want = want && (NXrow->Time-0.5*NXrow->TimeI >= sel->timeRange[0]);
    want = want && (NXrow->Time+0.5*NXrow->TimeI <= sel->timeRange[1]);

    /* If so write selected flagging */
    if (want) {
      nscan++; /* Count */
      FGrow->SourID       = NXrow->SourID;
      FGrow->SubA         = NXrow->SubA;
      FGrow->freqID       = NXrow->FreqID;
      FGrow->TimeRange[0] = NXrow->Time-0.51*NXrow->TimeI;
      FGrow->TimeRange[1] = NXrow->Time-0.5*NXrow->TimeI + begDrop;
      sumTime +=  begDrop;
     /* Beginning */
      if (begDrop> 0.0) {
	/* Loop over antennas */
	for (iant=0; iant<maxAnt; iant++) {
	  if (!AntOK[iant]) continue;  /* selected? */
	  /* write row */
	  FGrow->ants[0]      = iant+1;
	  orow = FGTab->myDesc->nrow+1;
	  ObitTableFGWriteRow (FGTab, orow, FGrow, err);
	  if (err->error) goto cleanup;
	}
      } /* end antenna loop */
      /* Loop over antennas */
	/* end */
      if (endDrop> 0.0) {
	FGrow->SourID       = NXrow->SourID;
	FGrow->SubA         = NXrow->SubA;
	FGrow->freqID       = NXrow->FreqID;
	FGrow->TimeRange[0] = NXrow->Time+0.5*NXrow->TimeI - endDrop;
	FGrow->TimeRange[1] = NXrow->Time+0.51*NXrow->TimeI;
	sumTime += endDrop;
	for (iant=0; iant<maxAnt; iant++) {
	  if (!AntOK[iant]) continue;  /* selected? */
	  /* write row */
	  FGrow->ants[0]      = iant+1;
	  orow = FGTab->myDesc->nrow+1;
	  ObitTableFGWriteRow (FGTab, orow, FGrow, err);
	  if (err->error) goto cleanup;
	}
      } /* end antenna loop */
    }
    
  } /* end loop over rows */

  /* close and cleanup */
 cleanup:
  if (AntOK) g_free (AntOK);
  ObitTableNXClose (in, err) ;
  ObitTableFGClose (FGTab, err) ;
  NXrow = ObitTableNXRowUnref(NXrow);
  FGrow = ObitTableFGRowUnref(FGrow);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Tell about it */
  Obit_log_error(err, OBIT_InfoErr, 
		 "Flagged %f min in %d scans", sumTime*1440.0, nscan);
  ObitErrLog(err); 

  return;
} /* end ObitTableNXUtilQuack */

