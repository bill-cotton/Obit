/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006,2016                                          */
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

#include "ObitTablePSUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTablePSUtil.c
 * ObitTablePS class utility function definitions.
 */

/*----------------------Public functions---------------------------*/


/**
 * Determine if a given source has successful processing indicated in the
 * Processing Summaty table.(STATUS=="Done" or "Failed")
 * \param in        PS table to access, opened and closed externally
 * \param target    Desired target field name.
 * \param theRow    [out] row number on which match found, even if not 
 *                  successful, else -1.
 * \param *err      ObitErr error stack.
 * \return FALSE if successful prior processing indicated, else TRUE
 */
gboolean ObitTablePSWantSour (ObitTablePS *in, gchar *target,
			      olong *theRow, ObitErr *err)
{
  gboolean out = TRUE;
  ObitTablePSRow *row = NULL;
  olong irow, mxlen;
  gchar *routine = "ObitTablePSWantSour";

  *theRow = -1;  /* If target not found */

  /* error checks */
  if (err->error) return out;

  row  = newObitTablePSRow (in);

  /* How many characters to check */
  mxlen = MIN (16, strlen(target));

  /* Loop over table checking for target */
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTablePSReadRow (in, irow, row, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
    if (row->status==-1) continue;  /* Deselected? */
    
    /* Is this it? */
    if (!strncmp (row->FieldName, target, mxlen)) {
      if (!strncmp (row->Status, "Done", 4))        out = FALSE;
      else if (!strncmp (row->Status, "Failed", 6)) out = TRUE;
      else out = TRUE;
      *theRow = irow;
      break;
    }
  } /* end loop over rows */

  /* release row object */
  row = ObitTablePSRowUnref(row);

  return out;
} /* end ObitTablePSWantSour */

/**
 * Write Processing Summary table from entries in an ObitInfoList
 * \param inUV      ObitUV object with PS table to update
 * \param info      List with parameters:
 * \li "FieldName" OBIT_string (16,1,1) 
 * \li "obsdate"   OBIT_string (8,1,1)  Observing ref. date.(YYYYMMDD)
 *                                      default - get from header
 * \li "RAPoint"   OBIT_double (1,1,1)  Pointing RA at epoch 2000 
 * \li "DecPoint"  OBIT_double (1,1,1)  Pointing Dec at epoch 2000
 * \li "Status"    OBIT_string (8,1,1)  Processing status, default "Done"
 * \li "IQU"       OBIT_bool   (3,1,1)  IQU Flags True if processed}
 * \li "Special"   OBIT_bool   (1,1,1)  Special processing applied, def. FALSE
 * \li "NoFields"  OBIT_int    (1,1,1)  Number of sub fields
 * \li "CCLimit"   OBIT_bool   (3,1,1)  Flags to using all possible CC, def FALSE
 * \li "PerCent"   OBIT_float  (3,1,1)  percent of visibility data I,Q,U, def 100.0
 * \li "PixMax"    OBIT_float  (3,1,1)  Pixel max (Jy) I,Q,U
 * \li "PixMin"    OBIT_float  (3,1,1)  Pixel min (Jy) I,Q,U
 * \li "Quality"   OBIT_float  (3,1,1)  Quality measure, RMS I,Q,U 
 * \li "Comment"   OBIT_string (80,1,1) Comments, default blank
 * \param *err      ObitErr error stack.
 */
void ObitTablePSSummary (ObitUV *inUV, ObitInfoList *info, ObitErr *err)
{
  ObitTablePS *PSTable=NULL;
  ObitTablePSRow *row=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong iver, theRow=-1;
  gchar FieldName[17], obsd[9], *cp;
  gchar *routine = "ObitTablePSSummary";

  /* error checks */
  if (err->error) return;

  /* Get/Create table */
  iver = 1;
  PSTable = newObitTablePSValue (inUV->name, (ObitData*)inUV, &iver, 
				 OBIT_IO_ReadWrite, err);
  if (PSTable) ObitTablePSOpen (PSTable, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* See if it's there and if so get row number */
  strcpy (FieldName, "Unspecified");
  ObitInfoListGetTest(info, "FieldName", &type, dim, FieldName);
  FieldName[16] = 0;
  ObitTablePSWantSour (PSTable, FieldName, &theRow, err);
 
  row  = newObitTablePSRow (PSTable);
  ObitTablePSSetRow (PSTable, row, err);

  /* If it exists, read old record */
  if (theRow>0) ObitTablePSReadRow (PSTable, theRow, row, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* Fill in values */
  strncpy (row->FieldName, FieldName, 16);
  /* Reorder obsdate to get rid of '-' */
  cp = inUV->myDesc->obsdat;
  if ((cp[4]=='-') && (cp[7]=='-')) {
    obsd[0] = cp[0]; obsd[1] = cp[1]; obsd[2] = cp[2]; obsd[3] = cp[3]; 
    obsd[4] = cp[5]; obsd[5] = cp[6]; obsd[6] = cp[8]; obsd[7] = cp[9]; 
  } else {
    strncpy (obsd, cp, 8);
  }
  strncpy (row->obsdate, obsd, 8);
  ObitInfoListGetTest(info, "obsdate", &type, dim, row->obsdate);
  strncpy (row->Status, "Done", 5);
  ObitInfoListGetTest(info, "Status", &type, dim, row->Status);
  strncpy (row->Comment, "                    ", 20);
  ObitInfoListGetTest(info, "Comment", &type, dim, row->Comment);
  row->RAPoint = 0.0;
  ObitInfoListGetTest(info, "RAPoint", &type, dim, &row->RAPoint);
  row->DecPoint = 0.0;
  ObitInfoListGetTest(info, "DecPoint", &type, dim, &row->DecPoint);
  row->IQU[0] = row->IQU[1] = row->IQU[2] = FALSE;
  ObitInfoListGetTest(info, "IQU", &type, dim, row->IQU);
  row->Special = FALSE;
  ObitInfoListGetTest(info, "Special", &type, dim, &row->Special);
  row->NoFields = 1;
  ObitInfoListGetTest(info, "NoFields", &type, dim, &row->NoFields);
  row->CCLimit[0] = row->CCLimit[1] = row->CCLimit[2] = FALSE;
  ObitInfoListGetTest(info, "CCLimit", &type, dim, row->CCLimit);
  row->PerCent[0] = row->PerCent[1] = row->PerCent[2] = 100.0;
  ObitInfoListGetTest(info, "PerCent", &type, dim, row->PerCent);
  row->PixMax[0] = row->PixMax[1] = row->PixMax[2] = 0.0;
  ObitInfoListGetTest(info, "PixMax", &type, dim, row->PixMax);
  row->PixMin[0] = row->PixMin[1] = row->PixMin[2] = 0.0;
  ObitInfoListGetTest(info, "PixMin", &type, dim, row->PixMin);
  row->Quality[0] = row->Quality[1] = row->Quality[2] = 0.0;
  ObitInfoListGetTest(info, "Quality", &type, dim, row->Quality);

  /* (re)Write row */
  ObitTablePSWriteRow (PSTable, theRow, row, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* Close up */
   ObitTablePSClose (PSTable,  err);
  if (PSTable) PSTable = ObitTablePSUnref(PSTable);   /* Done with table */

 /* release row object */
  row = ObitTablePSRowUnref(row);

} /* end ObitTablePSSummary */

