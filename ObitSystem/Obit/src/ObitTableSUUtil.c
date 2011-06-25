/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2010                                          */
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

#include "ObitTableSUUtil.h"
#include "ObitPrecess.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableSUUtil.c
 * ObitTableSU class utility function definitions.
 */

/*----------------------Public functions---------------------------*/

/**
 * Lookup a list of sources in an SU table.
 * This is intended for finding user selected sources.
 * \param in       Table to obtain data from
 * \param dim      dimensionality of inlist, first, the length if the source names, 
 *                 then the number of entries.
 * \param inlist   List of source names in single array with no nulls,
 *                 any leading '-' is stripped;
 *                 Any nonblank entries that are not found in the SU table
 *                 will generate  a warning.
 * \param Qual     Desired source qualifier, -1 => any
 * \param souCode  Source Cal code desired, '    ' => any code selected
 *                 '*   ' => any non blank code (calibrators only)
 *                 '-CAL' => blank codes only (no calibrators)
 * \param outlist  List of source IDs corresponding to inlist, -1 => not found.
 *                 Nothing done if this is NULL.
 * \param xselect  TRUE if source list selected, always TRUE as the list 
 *                 returned is those that match selection or don't match deselection.
 * \param Number   on input, the size allocated to outList, on output, the number found
 * \param *err     ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableSULookup (ObitTableSU *in, gint32 *dim, gchar *inlist, 
			      olong Qual, gchar souCode[5], olong *outlist, 
			      gboolean *xselect, olong *Number, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableSURow *row;
  olong i, j, l, maxNum, ncheck, *cross, size;
  gboolean select, gotSome, want, match;
  gchar tempName[101], temp2Name[101]; /* should always be big enough */
  gchar *routine = "ObitTableSULookup";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitTableSUIsA(in));
  if (outlist==NULL) return OBIT_IO_OK; /* Output list given */

  /* initialize */
  *xselect = TRUE;
  select   = TRUE;
  maxNum   = *Number;
  *Number  = 0;
  for (i=0; i<MIN (maxNum,dim[1]); i++) outlist[i] = -1;

  /* Open table */
  retCode = ObitTableSUOpen (in, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine,in->name, retCode);

  /* Create table row */
  row = newObitTableSURow (in);

  /* Cross list of sources */
  size  = MIN (maxNum,dim[1]);
  size  = MAX (size, in->myDesc->nrow);
  cross = g_malloc0((size+5)*sizeof(olong));
  for (i=0; i<size; i++) cross[i] = -1;

  /* loop over table */
  while (retCode==OBIT_IO_OK) {
    retCode = ObitTableSUReadRow (in, -1, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) {
      g_free(cross);
      Obit_traceback_val (err, routine,in->name, retCode);
    }
    /* Blank fill Source to dim[0], in temp2Name */
    for (i=0; i<dim[0]; i++) temp2Name[i] = ' ';  temp2Name[i] = 0;
    l = MIN (in->myDesc->repeat[in->SourceCol], strlen(row->Source));
    for (i=0; i<l; i++)  temp2Name[i] = row->Source[i];
    temp2Name[i] = 0;  /* to be sure terminated */

    /* loop through inlist and check if it matches */
    match = FALSE;
    for (i=0; i<dim[1]; i++) {
      for (j=0; j<dim[0]; j++) tempName[j] = ' '; tempName[j] = 0;
      /* get blank padded name from list */
      for (j=0; j<dim[0]; j++) {
	if (inlist[i*dim[0]+j]==0) break;  /* only values in string */
	tempName[j] = inlist[i*dim[0]+j]; 
      }
      /* blank fill to at least 16 characters */
      if (strlen(tempName)<16) {
	for (j=strlen(tempName); j<16; j++) tempName[j] = ' '; 
	tempName[j] = 0;
      }
      /* Have an initial '-'? */
      if (tempName[0]=='-') {
	select = FALSE;
	for (j=0; j<dim[0]; j++) tempName[j]= tempName[j+1]; /* remove '-' */
	tempName[dim[0]-1]=' ';tempName[dim[0]]=0;  /* Add blank to end */
      }
      /* Is this a match? */
      ncheck = MAX (strlen(tempName), strlen(temp2Name));
      match = (!strncmp (temp2Name, tempName, ncheck));
      if (match) break;
    }/* End loop over list */

    /* Match name and select true or, no match and select false, are OK */
    want = ((match) && (select)) || ((!match) && !(select));
    want = want && ((Qual==row->Qual) || (Qual<=-1));
    if (want && (strncmp (souCode, "    ", 4))) {
      if (!strncmp (souCode, "-CAL", 4))  /* want blank calcode */
	want = want && (!strncmp (row->CalCode, "    ", 4));
      else if (!strncmp (souCode, "*   ", 4)) /* want nonblank calcode */
	want = want && (strncmp (row->CalCode, "    ", 4));
      /* Must match */
      else want = want && (!strncmp (souCode, row->CalCode, 4));
    }
    if (want) {  /* Save number if wanted */
      if ((*Number)<maxNum) { /* array full? */
	outlist[(*Number)] = row->SourID;
	cross[i] = row->SourID;
	(*Number)++;   /* Count how many actually found */
      } else { /* Full */
	ObitTrimTrail(tempName);
	Obit_log_error(err, OBIT_InfoWarn, 
		       "%s: Warning: Source %s :%4.4d not added - table Full", 
		       routine, tempName, Qual);
      }
    } /* end if wanted */
  } /* End loop over table */
  
  /* check for errors */
  if ((retCode > OBIT_IO_EOF) || (err->error)) {
    g_free(cross);
    Obit_traceback_val (err, routine,in->name, retCode);
  }
  
  /* Release table row */
  row = ObitTableSURowUnref (row);

  /* Close table */
  retCode = ObitTableSUClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) {
    g_free(cross);
    Obit_traceback_val (err, routine,in->name, retCode);
  }

  /* Be sure all nonblank entries found else issue warning */
  gotSome = FALSE;  /* Any matches */
  for (i=0; i<dim[1]; i++) {
    gotSome = gotSome || (outlist[i]>=0);
    if (cross[i]<0) {/* Not found */
      /* get name from list */
      for (j=0; j<dim[0]; j++) tempName[j] = inlist[i*dim[0]+j]; tempName[j] = 0;
      /* Have an initial '-'? */
      if (tempName[0]=='-') for (j=0; j<dim[0]; j++) tempName[j]= tempName[j+1]; /* remove '-' */

      /* all blank is OK */
      if (strncmp(tempName, "                ", dim[0])) {
	ObitTrimTrail(tempName);
	Obit_log_error(err, OBIT_InfoWarn, 
		       "Source %s :%4.4d code %s not found in source table or dup in list", 
		       tempName, Qual, souCode);
      }
    }
  }

  /* Anything selected? */
  if (!gotSome) *xselect = FALSE;
  g_free(cross);  /* Free */

  return retCode;
} /* end ObitTableSULookup */

/**
 * Convert the contents of a ObitTableSU into an ObitSourceList.
 * \param in       Table to obtain data from
 * \param *err     ObitErr error stack.
 * \return requested ObitSourceList
 */
ObitSourceList* ObitTableSUGetList (ObitTableSU *in, ObitErr *err) {
  ObitSourceList *out=NULL;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableSURow *row;
  olong maxSUid, nsid, sid;
  olong i, irow;
  gchar tempName[101]; /* should always be big enough */
  gchar *routine = "ObitTableSUGetList";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitTableSUIsA(in));

  /* Open table */
  retCode = ObitTableSUOpen (in, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, out);
  
  /* Create table row */
  row = newObitTableSURow (in);

  /* loop over table looking for highest number */
  maxSUid = -1;
  nsid = 0;
  irow = 0;
  while (retCode==OBIT_IO_OK) {
    irow++;
    retCode = ObitTableSUReadRow (in, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, out);
    
    maxSUid = MAX (maxSUid, row->SourID);
    nsid++;
  } /* end loop over file */
  
  /* check for errors */
  if ((retCode > OBIT_IO_EOF) || (err->error))
    Obit_traceback_val (err, routine, in->name, out);

  /* Create output */
  g_snprintf (tempName, 100, "Source List for %s",in->name);
  out = ObitSourceListCreate (tempName, maxSUid);
  
  /* loop over table saving information */
  retCode = OBIT_IO_OK;
  irow = 0;
  while (retCode==OBIT_IO_OK) {
    irow++;
    retCode = ObitTableSUReadRow (in, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, out);

    /* Check for corruption by AIPS ++ */
    if (row->RAMean<0.0)   row->RAMean += 360.0;
    if (row->RAMean>360.0) row->RAMean -= 360.0;
    if (row->RAApp<0.0)    row->RAApp  += 360.0;
    if (row->RAApp>360.0)  row->RAApp  -= 360.0;

    sid = row->SourID - 1;
    out->SUlist[sid]->SourID  = row->SourID;
    out->SUlist[sid]->Qual    = row->Qual;
    out->SUlist[sid]->numIF   = in->numIF;
    out->SUlist[sid]->equinox = row->Epoch;   /* correct AIPS misnaming */
    out->SUlist[sid]->RAMean  = row->RAMean;
    /* Patch AIPS++ corruption */
    if (out->SUlist[sid]->RAMean<0.0) out->SUlist[sid]->RAMean += 360.0;
    out->SUlist[sid]->DecMean = row->DecMean;
    out->SUlist[sid]->RAApp   = row->RAApp;
    /* Patch AIPS++ corruption */
    if (out->SUlist[sid]->RAApp<0.0) out->SUlist[sid]->RAApp += 360.0;
    out->SUlist[sid]->DecApp  = row->DecApp;
    /* More Patch AIPS++ corruption */
    if ((out->SUlist[sid]->RAApp==0.0) && (out->SUlist[sid]->DecApp==0.0))
      ObitPrecessUVJPrecessApp (((ObitUV*)in->myHost)->myDesc, out->SUlist[sid]);
    out->SUlist[sid]->Bandwidth  = row->Bandwidth;
    strncpy (out->SUlist[sid]->SourceName, row->Source, 17);
    out->SUlist[sid]->SourceName[16] = 0;  /* to be sure */
    strncpy (out->SUlist[sid]->CalCode, row->CalCode, 5);
    out->SUlist[sid]->CalCode[4] = 0;  /* to be sure */
    if (out->SUlist[sid]->IFlux)     g_free (out->SUlist[sid]->IFlux);
    if (out->SUlist[sid]->QFlux)     g_free (out->SUlist[sid]->QFlux);
    if (out->SUlist[sid]->UFlux)     g_free (out->SUlist[sid]->UFlux);
    if (out->SUlist[sid]->VFlux)     g_free (out->SUlist[sid]->VFlux);
    if (out->SUlist[sid]->FreqOff)   g_free (out->SUlist[sid]->FreqOff);
    if (out->SUlist[sid]->LSRVel)    g_free (out->SUlist[sid]->LSRVel);
    if (out->SUlist[sid]->RestFreq)  g_free (out->SUlist[sid]->RestFreq);
    if (out->SUlist[sid]->numIF>0) {
      out->SUlist[sid]->IFlux     = g_malloc(out->SUlist[sid]->numIF*sizeof(ofloat));
      out->SUlist[sid]->QFlux     = g_malloc(out->SUlist[sid]->numIF*sizeof(ofloat));
      out->SUlist[sid]->UFlux     = g_malloc(out->SUlist[sid]->numIF*sizeof(ofloat));
      out->SUlist[sid]->VFlux     = g_malloc(out->SUlist[sid]->numIF*sizeof(ofloat));
      out->SUlist[sid]->FreqOff   = g_malloc(out->SUlist[sid]->numIF*sizeof(odouble));
      out->SUlist[sid]->LSRVel    = g_malloc(out->SUlist[sid]->numIF*sizeof(odouble));
      out->SUlist[sid]->RestFreq  = g_malloc(out->SUlist[sid]->numIF*sizeof(odouble));
      for (i=0; i<out->SUlist[sid]->numIF; i++)  out->SUlist[sid]->IFlux[i]    = row->IFlux[i];
      for (i=0; i<out->SUlist[sid]->numIF; i++)  out->SUlist[sid]->QFlux[i]    = row->QFlux[i];
      for (i=0; i<out->SUlist[sid]->numIF; i++)  out->SUlist[sid]->UFlux[i]    = row->UFlux[i];
      for (i=0; i<out->SUlist[sid]->numIF; i++)  out->SUlist[sid]->VFlux[i]    = row->VFlux[i];
      for (i=0; i<out->SUlist[sid]->numIF; i++)  out->SUlist[sid]->FreqOff[i]  = row->FreqOff[i];
      for (i=0; i<out->SUlist[sid]->numIF; i++)  out->SUlist[sid]->LSRVel[i]   = row->LSRVel[i];
      for (i=0; i<out->SUlist[sid]->numIF; i++)  out->SUlist[sid]->RestFreq[i] = row->RestFreq[i];
    }
    
    /* Save name as name of object */
    if (out->SUlist[sid]->name) g_free (out->SUlist[sid]->name);
    strncpy (tempName, row->Source, 100);
    tempName[in->myDesc->repeat[in->SourceCol]] = 0; /* Null terminate */
    out->SUlist[sid]->name = g_strdup(tempName);
    
  } /* end second loop over table */

  /* Dummy any missing sources */
  for (sid=0; sid<maxSUid; sid++) {
    if (out->SUlist[sid]->numIF<=0) {
      out->SUlist[sid]->numIF     = 0;
      out->SUlist[sid]->IFlux     = g_malloc0(sizeof(ofloat));
      out->SUlist[sid]->QFlux     = g_malloc0(sizeof(ofloat));
      out->SUlist[sid]->UFlux     = g_malloc0(sizeof(ofloat));
      out->SUlist[sid]->VFlux     = g_malloc0(sizeof(ofloat));
      out->SUlist[sid]->FreqOff   = g_malloc0(sizeof(odouble));
      out->SUlist[sid]->LSRVel    = g_malloc0(sizeof(odouble));
      out->SUlist[sid]->RestFreq  = g_malloc0(sizeof(odouble));
    }
  } /* end dummy missing Sources */

 /* Release table row */
  row = ObitTableSURowUnref (row);

  /* Close table */
  retCode = ObitTableSUClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, out);

  return out;
} /* end ObitTableSUGetList */

/**
 * Copies SU tables from inUV to outUV with selection in inUV
 * Copies all sources but selects IFs
 * \param inUV     Input UV to copy from
 * \param outUV    Output UV to copy to
 * \param err      ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableSUSelect (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableSU    *inTab=NULL, *outTab=NULL;
  ObitTableSURow *inRow=NULL, *outRow=NULL;
  ObitSource *source=NULL;
  olong i, iif, oif;
  olong highSUver, iSUver, inSURow, outSURow;
  oint numIF;
  gchar *SUType = "AIPS SU";
  gchar *routine = "ObitTableSUSelect";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));

  /* Fully instantiate UV files */
  ObitUVFullInstantiate (inUV, TRUE, err);
  if (err->error )Obit_traceback_val (err, routine, inUV->name, retCode);
  ObitUVFullInstantiate (outUV, FALSE, err);
  if (err->error )Obit_traceback_val (err, routine, outUV->name, retCode);

  /* How many SU tables  */
  highSUver = ObitTableListGetHigh (inUV->tableList, SUType);

  /* Are there any? */
  if (highSUver <= 0) return OBIT_IO_OK;

  /* Only 1 SU table allowed */
  iSUver = 1;

  /* Get input table */
  numIF = 0;
  inTab = 
    newObitTableSUValue (inUV->name, (ObitData*)inUV, &iSUver, OBIT_IO_ReadOnly, 
			 numIF, err);
  if (err->error) Obit_traceback_val (err, routine, inTab->name, retCode);
  /* Find it */
  if (inTab==NULL) return OBIT_IO_OK; /* No - return */
  
  /* Open input table */
  retCode = ObitTableSUOpen (inTab, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inTab->name, retCode);
  
  /* Delete any old output table */
  retCode = ObitDataZapTable ((ObitData*)outUV, SUType, iSUver, err);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, retCode);
  
  /* Create output table */
  numIF  = inUV->mySel->numberIF;
  outTab = 
    newObitTableSUValue (outUV->name, (ObitData*)outUV, &iSUver, OBIT_IO_WriteOnly, 
			 numIF, err);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, retCode);
  /* Create it? */
  Obit_retval_if_fail((outTab!=NULL), err, retCode,
		      "%s: Could not create SU table %d for %s", 
		      routine, iSUver, outTab->name);
  
  /* Open output table */
  retCode = ObitTableSUOpen (outTab, OBIT_IO_WriteOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inTab->name, retCode);
  
  /* Update header info */
  for (i=0; i<MAXKEYCHARTABLESU; i++)  outTab->velType[i] = inTab->velType[i];
  for (i=0; i<MAXKEYCHARTABLESU; i++)  outTab->velDef[i]  = inTab->velDef[i];
  outTab->FreqID    = inTab->FreqID;
  
  /* Set rows */
  inRow  = newObitTableSURow (inTab);
  outRow = newObitTableSURow (outTab);
  ObitTableSUSetRow (outTab, outRow, err);
  if (err->error) Obit_traceback_val (err, routine, outTab->name, retCode);
  
  /* Loop over table copying selected data */
  outSURow = -1;
  for (inSURow=1; inSURow<=inTab->myDesc->nrow; inSURow++) {
    retCode = ObitTableSUReadRow (inTab, inSURow, inRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inUV->name, retCode);
    if (inRow->status==-1) continue;
    
    /* Copy selected data */
    outRow->SourID    = inRow->SourID;
    outRow->Qual      = inRow->Qual;
    outRow->Bandwidth = inRow->Bandwidth;
    outRow->RAMean    = inRow->RAMean;
    if (outRow->RAMean<0.0) outRow->RAMean += 360.0; /* Patch AIPS++ corruption */
    outRow->DecMean   = inRow->DecMean;
    outRow->Epoch     = inRow->Epoch;
    outRow->RAApp     = inRow->RAApp;
    if (outRow->RAApp<0.0) outRow->RAApp += 360.0; /* Patch AIPS++ corruption */
    outRow->DecApp    = inRow->DecApp;
    /* Patch AIPS++ bug */
    if ((outRow->RAApp==0.0) && (outRow->DecApp==0.0)) {
      source = newObitSource("Temp");
      source->equinox = outUV->myDesc->equinox;
      source->RAMean  = outUV->myDesc->crval[outUV->myDesc->jlocr];
      source->DecMean = outUV->myDesc->crval[outUV->myDesc->jlocd];
      /* Compute apparent position */
      ObitPrecessUVJPrecessApp (outUV->myDesc, source);
      outRow->RAApp  = source->RAApp;
      outRow->DecApp = source->DecApp;
      source = ObitSourceUnref(source);
    }
    outRow->PMRa      = inRow->PMDec;
    for (i=0; i<inTab->myDesc->repeat[inTab->SourceCol]; i++) 
      outRow->Source[i] = inRow->Source[i];
    for (i=0; i<inTab->myDesc->repeat[inTab->CalCodeCol]; i++) 
      outRow->CalCode[i] = inRow->CalCode[i];

    oif = 0;
    for (iif=inUV->mySel->startIF-1; 
	 iif<inUV->mySel->startIF+inUV->mySel->numberIF-1;
	 iif++) {
      outRow->IFlux[oif]    = inRow->IFlux[iif];
      outRow->QFlux[oif]    = inRow->QFlux[iif];
      outRow->UFlux[oif]    = inRow->UFlux[iif];
      outRow->VFlux[oif]    = inRow->VFlux[iif];
      outRow->FreqOff[oif]  = inRow->FreqOff[iif];
      outRow->LSRVel[oif]   = inRow->LSRVel[iif];
      outRow->RestFreq[oif] = inRow->RestFreq[iif];
      oif++;
    }
    retCode = ObitTableSUWriteRow (outTab, outSURow, outRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inUV->name, retCode);
  } /* end loop over rows */
  
  /* Close tables */
  retCode = ObitTableSUClose (inTab, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inTab->name, retCode);
  retCode = ObitTableSUClose (outTab, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, outTab->name, retCode);

  /* release table objects */
  inTab  = ObitTableSUUnref(inTab);
  outTab = ObitTableSUUnref(outTab);
  
  /* release row objects */
  inRow  = ObitTableSURowUnref(inRow);
  outRow = ObitTableSURowUnref(outRow);

  return retCode;
} /* end ObitTableSUSelect */

