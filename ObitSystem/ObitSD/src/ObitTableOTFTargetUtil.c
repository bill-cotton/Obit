/* $Id: ObitTableOTFTargetUtil.c,v 1.7 2008/02/28 15:23:46 bcotton Exp $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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

#include "ObitTableOTFTargetUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableOTFTargetUtil.c
 * ObitTableOTFTarget class utility function definitions.
 */

/*----------------------Public functions---------------------------*/

/**
 * Lookup a list of sources in an OTFTarget table.
 * This is intended for finding user selected sources.
 * \param in       Table to obtain data from
 * \param dim      dimensionality of inlist, first, the length of the source names, 
 *                 then the number of entries.
 * \param inlist   List of source names in single array with no nulls,
 *                 Any nonblank entries that are not found in the OTFTarget table
 *                 will generate  a warning.
 * \param outlist  List of source IDs corresponding to inlist, -1 => not found.
 * \param *err     ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableOTFTargetLookup (ObitTableOTFTarget *in, gint32 *dim, 
				     gchar *inlist, olong *outlist, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableOTFTargetRow *row;
  olong i, j;
  gchar tempName[101],tempName2[101] ; /* should always be big enough */
  gchar *routine = "ObitTableOTFTargetGetLookup";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitTableOTFTargetIsA(in));

  /* initialize */
  for (i=0; i<dim[1]; i++) outlist[i] = -1;

  /* Open table */
  retCode = ObitTableOTFTargetOpen (in, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Create table row */
  row = newObitTableOTFTargetRow (in);

  /* loop over table */
  while (retCode==OBIT_IO_OK) {
    retCode = ObitTableOTFTargetReadRow (in, -1, row, err);
    if (retCode == OBIT_IO_EOF) break; 
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, retCode);

    /* loop through inlist and check if it matches */
    for (i=0; i<dim[1]; i++) {
      /* get name from list */
      for (j=0; j<dim[0]; j++) tempName[j] = inlist[i*dim[0]+j]; tempName[j] = 0;
      /* drop trailing blanks */
      for (j=strlen(tempName)-1; j>=0; j--) {
	if (tempName[j]==' ') tempName[j] = 0;
	else break;
      } 
      /* Name from table */
      strncpy (tempName2, row->Target, 100);
       /* drop trailing blanks */
      for (j=in->myDesc->repeat[in->TargetCol]-1; j>=0; j--) {
	if (tempName2[j]==' ') tempName2[j] = 0;
	else break;
      } 
     /* Is this a match? */
      if ((!strncmp (tempName2, tempName, dim[0])) && 
	  (strlen(tempName)==strlen(tempName2)))
	outlist[i] = row->TargID;
    }
  }
  
  /* check for errors */
  if ((retCode > OBIT_IO_EOF) || (err->error))
    Obit_traceback_val (err, routine,in->name, retCode);
  
  /* Release table row */
  row = ObitTableOTFTargetRowUnref (row);

  /* Close table */
  retCode = ObitTableOTFTargetClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine,in->name, retCode);

  /* Be sure all nonblank entries found else issue warning */
  for (i=0; i<dim[1]; i++) {
    if (outlist[i]<0) {/* Not found */
      /* get name from list */
      for (j=0; j<dim[0]; j++) tempName[j] = inlist[i*dim[0]+j]; tempName[j] = 0;
      /* all blank is OK */
      if (strncmp(tempName, "                ", dim[0]))
	  Obit_log_error(err, OBIT_InfoWarn, 
			 "Warning: Target %s not found in target table", tempName);
    }
  }
  return retCode;
} /* end ObitTableOTFTargetLookup */

/**
 * Convert the contents of a ObitTableOTFTarget into an ObitSourceList.
 * \param in       Table to obtain data from
 * \param *err     ObitErr error stack.
 * \return requested ObitSourceList
 */
ObitSourceList* ObitTableOTFTargetGetList (ObitTableOTFTarget *in, ObitErr *err) {
  ObitSourceList *out=NULL;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableOTFTargetRow *row;
  olong maxTargetID, sid;
  olong irow;
  gchar tempName[101]; /* should always be big enough */
  gchar *routine = "ObitTableOTFTargetGetList";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitTableOTFTargetIsA(in));

  /* Open table */
  retCode = ObitTableOTFTargetOpen (in, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, out);
  
  /* Create table row */
  row = newObitTableOTFTargetRow (in);

  /* loop over table looking for highest number */
  maxTargetID = -1;
  irow = 0;
  while (retCode==OBIT_IO_OK) {
    irow++;
    retCode = ObitTableOTFTargetReadRow (in, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, out);
    
    maxTargetID = MAX (maxTargetID, row->TargID);

  } /* end loop over file */
  
  /* check for errors */
  if ((retCode > OBIT_IO_EOF) || (err->error))
    Obit_traceback_val (err, routine, in->name, out);

  /* Create output */
  g_snprintf (tempName, 100, "Target List for %s",in->name);
  out = ObitSourceListCreate (tempName, maxTargetID);
  
  /* loop over table saving information */
  retCode = OBIT_IO_OK;
  irow = 0;
  while (retCode==OBIT_IO_OK) {
    irow++;
    retCode = ObitTableOTFTargetReadRow (in, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, out);

    sid = row->TargID - 1;
    out->SUlist[sid]->SourID  = row->TargID;
    out->SUlist[sid]->RAMean  = row->RAMean;
    out->SUlist[sid]->Qual    = row->Qual;
    out->SUlist[sid]->DecMean = row->DecMean;
    out->SUlist[sid]->RAApp   = row->RAApp;
    out->SUlist[sid]->DecApp  = row->DecApp;

    /* Save name as name of object */
    if (out->SUlist[sid]->name) g_free (out->SUlist[sid]->name);
    strncpy (tempName, row->Target, 100);
    tempName[in->myDesc->repeat[in->TargetCol]] = 0; /* Null terminate */
    out->SUlist[sid]->name = g_strdup(tempName);
    strncpy (out->SUlist[sid]->SourceName, tempName, 20);

  } /* end second loop over table */

 /* Release table row */
  row = ObitTableOTFTargetRowUnref (row);

  /* Close table */
  retCode = ObitTableOTFTargetClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, out);

  return out;
} /* end ObitTableOTFTargetGetList */

/**
 * Get source info for a given source ID
 * \param in       Table to obtain data from
 * \param targID   target ID to look up
 * \param RA       [Out] RA of mean epoch
 * \param Dec      [Out] Dec of mean epoch
 * \param Flux     [Out] Flux density
 * \param *err     ObitErr error stack.
 * \return requested ObitSourceList
 */
ObitIOCode ObitTableOTFTargetGetSource (ObitTableOTFTarget *in, olong targID, 
					odouble *RA, odouble *Dec, ofloat *Flux, 
					ObitErr *err) 
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableOTFTargetRow *row;
  olong irow;
  gchar *routine = "ObitTableOTFTargetGetList";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitTableOTFTargetIsA(in));

  /* Open table */
  retCode = ObitTableOTFTargetOpen (in, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, retCode);
  
  /* Create table row */
  row = newObitTableOTFTargetRow (in);

  /* loop over table looking for targID */
  irow = 0;
  while (retCode==OBIT_IO_OK) {
    irow++;
    retCode = ObitTableOTFTargetReadRow (in, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, retCode);
    if (row->TargID==targID) break;  /* found it?*/
  } /* end loop over file */

  /* Found it? */
  if (row->TargID!=targID) {
    Obit_log_error(err, OBIT_InfoErr,
		   "Target %d not found in target table", targID);
    return OBIT_IO_SpecErr;
  }
  
  /* Save info */
  *RA  = row->RAMean;
  *Dec = row->DecMean;
  *Flux = row->IFlux;

 /* Release table row */
  row = ObitTableOTFTargetRowUnref (row);

  /* Close table */
  retCode = ObitTableOTFTargetClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, retCode);

  return retCode;
} /* end ObitTableOTFTargetGetSource */

/**
 * Set source info for a given source ID
 * \param in       Table to update
 * \param targID   target ID to look up
 * \param RA       RA of mean epoch
 * \param Dec      Dec of mean epoch
 * \param Flux     Flux density
 * \param *err     ObitErr error stack.
 * \return requested ObitSourceList
 */
ObitIOCode ObitTableOTFTargetSetSource (ObitTableOTFTarget *in, olong targID, 
					odouble RA, odouble Dec, ofloat Flux, 
					ObitErr *err) 
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableOTFTargetRow *row;
  olong irow;
  gchar *routine = "ObitTableOTFTargetSetList";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitTableOTFTargetIsA(in));

  /* Open table */
  retCode = ObitTableOTFTargetOpen (in, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, retCode);
  
  /* Create table row */
  row = newObitTableOTFTargetRow (in);

  /* loop over table looking for targID */
  irow = 0;
  while (retCode==OBIT_IO_OK) {
    irow++;
    retCode = ObitTableOTFTargetReadRow (in, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, retCode);
    if (row->TargID==targID) break;  /* found it?*/
  } /* end loop over file */

  /* Found it? */
  if (row->TargID!=targID) {
    Obit_log_error(err, OBIT_InfoErr,
		   "Target %d not found in target table", targID);
    return OBIT_IO_SpecErr;
  }
  
  /* Save info */
  row->RAMean  = RA;
  row->DecMean = Dec;
  row->IFlux   = Flux;

  /* rewrite row */
  retCode = ObitTableOTFTargetWriteRow (in, irow, row, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Close table */
  retCode = ObitTableOTFTargetClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, retCode);

 /* Release table row */
  row = ObitTableOTFTargetRowUnref (row);

  return retCode;
} /* end ObitTableOTFTargetSetSource */

/**
 * Determine target number, adding new entry if necessary
 * \param in       Table to update
 * \param name     target name to look up, should be NULL terminated
 * \param qual     target qualifier to look up
 * \param RA       RA of mean epoch
 * \param Dec      Dec of mean epoch
 * \param equinox  Equinox of the position
 * \param *err     ObitErr error stack.
 * \return requested target ID
 */
olong ObitTableOTFTargetGetAddSource (ObitTableOTFTarget *in, gchar *name, olong qual,
				     odouble RA, odouble Dec, odouble equinox, 
				     ObitErr *err)
{
  olong targ=0;
  ObitTableOTFTargetRow* row;
  olong iRow;
  gboolean isNew,doWrite;
  gchar tName[24], t2Name[24];
  gchar *routine = "ObitTableOTFTargetGetAddSource";

  strncpy (t2Name, name, 23); t2Name[23]=0;
  ObitTrimTrail(t2Name);

  /* Open table */
  if ((ObitTableOTFTargetOpen (in, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output OTFTarget table");
    return targ;
  }

  /* Create Row */
  row = newObitTableOTFTargetRow (in);

  /* attach to table buffer */
  ObitTableOTFTargetSetRow (in, row, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, targ);

  /* Newly created?  Just write new one */
  doWrite = FALSE;
  isNew   = (in->myDesc->nrow<=0);
  if (isNew) {
    targ = 1;
    row->TargID = targ;
    row->Qual   = qual;
    strncpy(row->Target, t2Name, 16);
    row->RAMean  = RA;
    row->DecMean = Dec;
    row->Epoch   = equinox;
    doWrite      = TRUE;
  } else { /* Existing, see if already exists? */

    /* loop through table */
    for (iRow = 1; iRow<=in->myDesc->nrow; iRow++) {
      if ((ObitTableOTFTargetReadRow (in, iRow, row, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading OTFTarget Table file");
	return targ;
      }
      strncpy (tName, row->Target, 16); tName[16]=0;
      ObitTrimTrail(tName);
      if ((!strncmp (tName, t2Name, 16)) && (row->Qual==qual)) {
	/* Found match */
	targ = row->TargID;
	break;
      }  
    } /* end loop over table */

    /* Add new entry? */
    if (targ<=0) {
      targ = in->myDesc->nrow + 1;
      row->TargID  = targ;
      row->Qual    = qual;
      row->RAMean  = RA;
      row->DecMean = Dec;
      row->Epoch   = equinox;
      strncpy(row->Target, name, 16);
      doWrite = TRUE;
    }
  } /* end output table already exists */

  /* need to write new entry? */
  if (doWrite) {
    iRow = in->myDesc->nrow + 1;
    if ((ObitTableOTFTargetWriteRow (in, iRow, row, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR writing OTFTarget Table file");
      return targ;
    }
  }
  
 /* Close  table */
  if ((ObitTableOTFTargetClose (in, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output OTFTarget Table file");
    return targ;
  }

  /* Cleanup */
  row = ObitTableOTFTargetRowUnref(row);

  return targ;
} /* end ObitTableOTFTargetGetAddSource */

/**
 * Determine target number, position from name
 * \param in       Table to search
 * \param name     target name to look up, should be NULL terminated
 * \param qual     target qualifier to look up, 0=> any
 * \param RA       [out] RA of mean epoch
 * \param Dec      [out] Dec of mean epoch
 * \param *err     ObitErr error stack.
 * \return requested target ID
 */
olong ObitTableOTFTargetGetByName (ObitTableOTFTarget *in, gchar *name, olong qual,
				  odouble *RA, odouble *Dec, ObitErr *err)
{
  olong targ=0;
  ObitTableOTFTargetRow* row;
  olong iRow;
  gchar tName[24], t2Name[24];
  /*gchar *routine = "ObitTableOTFTargetGetByName";*/

  strncpy (t2Name, name, 23); t2Name[23]=0;
  ObitTrimTrail(t2Name);

  /* Open table */
  if ((ObitTableOTFTargetOpen (in, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output OTFTarget table");
    return targ;
  }

  /* Create Row */
  row = newObitTableOTFTargetRow (in);
  
  /* loop through table */
  for (iRow = 1; iRow<=in->myDesc->nrow; iRow++) {
    if ((ObitTableOTFTargetReadRow (in, iRow, row, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR reading OTFTarget Table file");
      return targ;
    }
    strncpy (tName, row->Target, 16); tName[16]=0;
    ObitTrimTrail(tName);
    if ((!strncmp (tName, t2Name, 16)) && (row->Qual==qual)) {
      /* Found match */
      targ = row->TargID;
      *RA  = row->RAMean;
      *Dec = row->DecMean;
      break;
    }  
  } /* end loop over table */
  
  /* Close  table */
  if ((ObitTableOTFTargetClose (in, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output OTFTarget Table file");
    return targ;
  }
  
  /* Cleanup */
  row = ObitTableOTFTargetRowUnref(row);
  
  return targ;
} /* end ObitTableOTFTargetGetByName */
