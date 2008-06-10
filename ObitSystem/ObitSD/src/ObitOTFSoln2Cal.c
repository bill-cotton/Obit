/* $Id: ObitOTFSoln2Cal.c,v 1.7 2008/02/13 21:13:13 bcotton Exp $ */
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

#include "ObitOTFSoln2Cal.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFSoln2Cal.c
 * ObitOTFSoln2Cal module calibration function definitions.
 */

/*---------------Private function prototypes----------------*/
/* Copy an open Soln table to an open Cal table */
void ObitOTFSolnCopyCal (ObitTableOTFSoln *inSoln, ObitTableOTFCal *outCal, 
			 ObitErr *err);

/*----------------------Public functions---------------------------*/
/**
 * Apply a Soln table to a Cal table and write a new Cal table
 * Calibration parameters are on the inOTF info member.
 * If an input Cal table is given then apply Solutions in this routine,
 * if no input Cal table, then copy the Soln table to a new Cal table
 * in ObitOTFSolnCopyCal.
 * \li "solnUse"   OBIT_int (1,1,1) Input Solution table version 
 * \li "calIn"     OBIT_int (1,1,1) Input Cal table version 
 *                 iff <0 then no input Cal table, copy Soln records to output.
 * \li "calOut"    OBIT_int (1,1,1) Output Calibration table version
 * \param inOTF    Input OTF data. 
 * \param outOTF   OTF with which the output  OTFCal is to be associated
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created OTFCal object which is 
 * associated with inOTF.
 */
ObitTableOTFCal* ObitOTFSoln2Cal (ObitOTF *inOTF, ObitOTF *outOTF, ObitErr *err)
{
  ObitTableOTFSoln *inSoln=NULL;
  ObitTableOTFSolnRow *inSolnRow=NULL;
  ObitTableOTFCal *outCal=NULL, *inCal=NULL;
  ObitTableOTFCalRow *outCalRow=NULL, *inCalRow=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat Wtp, Wtf, minMult, fblank = ObitMagicF();
  olong iRow, sRow, ver, i, j, npoly, inpoly, ndetect, limitS, limitC;
  olong  solnuse, calin, calout, iPr, iFl;
  ObitIOCode retCode;
  ObitOTFArrayGeom* geom;
  ofloat SolnTime[2], SolndAz[2], SolndEl[2];
  ofloat *SolnCal[2], *SolnAdd[2], *SolnMult[2], *SolnWt[2], *SolnPoly[2];
  gchar *tname;
  gchar *routine = "ObitOTFSoln2Cal";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outCal;
  g_assert (ObitOTFIsA(inOTF));
  g_assert (ObitOTFArrayGeomIsA(inOTF->geom));

  /* open/close OTF data to fully instantiate if not already open */
  if ((inOTF->myStatus==OBIT_Inactive) || (inOTF->myStatus==OBIT_Defined)) {
    retCode = ObitOTFOpen (inOTF, OBIT_IO_ReadWrite, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);
  }
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);


  /* Get array geometry */
  geom = inOTF->geom;

  /* Get table numbers */
  solnuse = 0;
  ObitInfoListGetTest(inOTF->info, "solnUse", &type, dim, (gpointer*)&solnuse);
  calin = 0;
  ObitInfoListGetTest(inOTF->info, "calIn", &type, dim, (gpointer*)&calin);
  calout = 0;
  ObitInfoListGetTest(inOTF->info, "calOut", &type, dim, (gpointer*)&calout);

  /* Instantiate input solution table object */
  tname = g_strconcat ("Input Solution for: ",inOTF->name, NULL);
  ver = solnuse;
  inSoln = newObitTableOTFSolnValue(tname, (ObitData*)inOTF, &ver, OBIT_IO_ReadWrite,  
				   inOTF->geom->numberDetect, 0, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);
 
  /* Open table */
  if ((ObitTableOTFSolnOpen (inSoln, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening input OTFSoln table");
    return outCal;
  }
  inSoln->numPoly = inSoln->myDesc->repeat[inSoln->polyCol];
  inSoln->numDet  = inSoln->myDesc->repeat[inSoln->addCol];

  /* Check that there are entries */
  if (inSoln->myDesc->nrow<=0) {
    Obit_log_error(err, OBIT_Error, "Input OTFSoln table has no data");
    return outCal;
  }

  /* if input Cal table exists then open */
  if (calin>=0) {
    /* Instantiate input calibration table object if one exists */
    tname = g_strconcat ("Input Calibration for: ",inOTF->name, NULL);
    ver = calin;
    inCal = newObitTableOTFCalValue(tname, (ObitData*)outOTF, &ver, OBIT_IO_ReadWrite,  
				    0, 0, err);
    g_free (tname);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);
    
    /* Open table */
    if ((ObitTableOTFCalOpen (inCal, OBIT_IO_ReadWrite, err) 
	 != OBIT_IO_OK) || (err->error))  { /* error test */
      Obit_log_error(err, OBIT_Error, "ERROR opening input OTFCal table");
      return outCal;
    }
    
    /* Better have some data */
    if (inCal->myDesc->nrow<=0) {
      Obit_log_error(err, OBIT_Error, "Input OTFCal table has no data");
      return outCal;
    }
  } /* end open input OTFCal table */

  /* Instantiate/Create output */
  tname = g_strconcat ("Output Calibration for: ",inOTF->name, NULL);
  ver = calout;
  if (calin>=0) { /* input Cal table? */
    inpoly  = inCal->myDesc->repeat[inCal->polyCol]; /* number of polynomial terms */
  } else { /* Nope - only Soln */
     inpoly  = inSoln->myDesc->repeat[inSoln->polyCol]; /* number of polynomial terms */
  }
  outCal = newObitTableOTFCalValue(tname, (ObitData*)outOTF, &ver, OBIT_IO_WriteOnly,  
				   inOTF->geom->numberDetect, inpoly, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);

  /* Open table */
  if ((ObitTableOTFCalOpen (outCal, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output Calibration table");
    return outCal;
  }

  /* if no input Cal table exists then copy input Soln table to output Cal table */
  if (calin<0) {
    ObitOTFSolnCopyCal (inSoln, outCal, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);
   return outCal;  /* done */
 }

  /* Create Soln Row */
  inSolnRow = newObitTableOTFSolnRow (inSoln);

  /* Create input Cal Row */
  inCalRow = newObitTableOTFCalRow (inCal);
  
  /* Attach input Cal row to output buffer */
  ObitTableOTFCalSetRow (inCal, inCalRow, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);
  
  /* Create output Row */
  outCalRow = newObitTableOTFCalRow (outCal);

  /* Attach row to output buffer */
  ObitTableOTFCalSetRow (outCal, outCalRow, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);

  /* how many of things? */
  inpoly  = inCal->myDesc->repeat[inCal->polyCol];
  npoly   = inSoln->myDesc->repeat[inSoln->polyCol];
  ndetect = inSoln->myDesc->repeat[inSoln->addCol];

  /* Create arrays to store Soln solutions to be interpolated */
  SolnCal[0]  = g_malloc0(ndetect*sizeof(ofloat));
  SolnAdd[0]  = g_malloc0(ndetect*sizeof(ofloat));
  SolnMult[0] = g_malloc0(ndetect*sizeof(ofloat));
  SolnWt[0]   = g_malloc0(ndetect*sizeof(ofloat));
  SolnPoly[0] = g_malloc0(inpoly*sizeof(ofloat));
  SolnCal[1]  = g_malloc0(ndetect*sizeof(ofloat));
  SolnAdd[1]  = g_malloc0(ndetect*sizeof(ofloat));
  SolnMult[1] = g_malloc0(ndetect*sizeof(ofloat));
  SolnPoly[1] = g_malloc0(inpoly*sizeof(ofloat));
  SolnWt[1]   = g_malloc0(ndetect*sizeof(ofloat));
 
  /* Initialize output row */
  outCalRow->dAz = 0.0;  /* no pointing corrections */
  outCalRow->dEl = 0.0; /* no pointing corrections */
  outCalRow->Target = 0;
  for (j=0; j<npoly; j++)   outCalRow->poly[j] = 0.0;
  for (j=0; j<ndetect; j++) outCalRow->mult[j] = 1.0;
  for (j=0; j<ndetect; j++) outCalRow->wt[j]   = 1.0;
  for (j=0; j<ndetect; j++) outCalRow->add[j]  = 0.0;
  for (j=0; j<ndetect; j++) outCalRow->cal[j]  = 0.0;

  retCode = OBIT_IO_OK;
  limitS = inSoln->myDesc->nrow; /* how many Soln rows to read */
  limitC = inCal->myDesc->nrow;  /* how many Cal rows to read */

  /* read first two Soln rows */
  sRow = 1;
  iPr = 0;  /* prior solution index */
  retCode = ObitTableOTFSolnReadRow (inSoln, sRow, inSolnRow, err);
  if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);
  /* Copy to storage */
  SolnTime[iPr] = inSolnRow->Time;
  SolndAz[iPr]  = inSolnRow->dAz;
  SolndEl[iPr] = inSolnRow->dEl;
  for (j=0; j<ndetect; j++)  SolnCal[iPr][j]  = inSolnRow->cal[j];
  for (j=0; j<ndetect; j++)  SolnAdd[iPr][j]  = inSolnRow->add[j];
  for (j=0; j<ndetect; j++)  SolnMult[iPr][j] = inSolnRow->mult[j];
  for (j=0; j<ndetect; j++)  SolnWt[iPr][j]   = inSolnRow->wt[j];
  for (j=0; j<npoly; j++)    SolnPoly[iPr][j] = inSolnRow->poly[j];

  if (limitS>1) {
    sRow++;
    iFl = 1;  /* following solution index */
    retCode = ObitTableOTFSolnReadRow (inSoln, sRow, inSolnRow, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);
    /* Copy to storage */
    SolnTime[iFl] = inSolnRow->Time;
    SolndAz[iFl]  = inSolnRow->dAz;
    SolndEl[iFl] = inSolnRow->dEl;
    for (j=0; j<ndetect; j++)  SolnCal[iFl][j]  = inSolnRow->cal[j];
    for (j=0; j<ndetect; j++)  SolnAdd[iFl][j]  = inSolnRow->add[j];
    for (j=0; j<ndetect; j++)  SolnMult[iFl][j] = inSolnRow->mult[j];
    for (j=0; j<ndetect; j++)  SolnWt[iFl][j]   = inSolnRow->wt[j];
    for (j=0; j<npoly; j++)    SolnPoly[iFl][j] = inSolnRow->poly[j];

  } else { /* only one Soln entry */
    iFl = iPr;
  }

  /* Loop interpolating Soln table and applying to input Cal table */
  Wtp = 1.0;
  Wtf = 0.0;
  for (iRow=0; iRow<=limitC; iRow++) {

    /* read input row */
    retCode = ObitTableOTFCalReadRow (inCal, iRow, inCalRow, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);
    if (retCode==OBIT_IO_EOF) break; /* done? */

    /* Copy input to output row */
    outCalRow->Time   = inCalRow->Time; 
    outCalRow->TimeI  = inCalRow->TimeI; 
    outCalRow->Target = inCalRow->Target; 
    outCalRow->dAz    = inCalRow->dAz; 
    outCalRow->dEl   = inCalRow->dEl;
    outCalRow->Target = inCalRow->Target;
    for (j=0; j<inpoly; j++)  outCalRow->poly[j] = inCalRow->poly[j];
    for (j=0; j<ndetect; j++) outCalRow->mult[j] = inCalRow->mult[j];
    for (j=0; j<ndetect; j++) outCalRow->wt[j]   = inCalRow->wt[j];
    for (j=0; j<ndetect; j++) outCalRow->add[j]  = inCalRow->add[j];
    for (j=0; j<ndetect; j++) outCalRow->cal[j]  = inCalRow->cal[j];

    /* Find solution bracketing the Cal entry */
    while ((outCalRow->Time > SolnTime[iFl]) && (sRow<=limitS)) {
      /* Need to read more solutions */
      sRow++;
      iPr = iFl;
      iFl = 1 - iFl;  /* swap solutions  */
      retCode = ObitTableOTFSolnReadRow (inSoln, sRow, inSolnRow, err);
      if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);
      /* Copy to storage */
      SolnTime[iFl] = inSolnRow->Time;
      SolndAz[iFl] = inSolnRow->dAz;
      SolndEl[iFl] = inSolnRow->dEl;
      for (j=0; j<ndetect; j++)  SolnCal[iFl][j]  = inSolnRow->cal[j];
      for (j=0; j<ndetect; j++)  SolnAdd[iFl][j]  = inSolnRow->add[j];
      for (j=0; j<ndetect; j++)  SolnMult[iFl][j] = inSolnRow->mult[j];
      for (j=0; j<ndetect; j++)  SolnWt[iFl][j  ] = inSolnRow->wt[j];
      for (j=0; j<npoly; j++)    SolnPoly[iFl][j] = inSolnRow->poly[j];
    }

    /* Set weights for linear interpolation */
    /* Cal entry before the first Soln entry? */
    if (outCalRow->Time < SolnTime[iPr]) {
      Wtp = 1.0;
      Wtf = 0.0;

      /* Between the two soln entries? */
    } else if ((outCalRow->Time >= SolnTime[iPr]) && 
	       (outCalRow->Time <= SolnTime[iFl])) {
      Wtf = (outCalRow->Time - SolnTime[iPr]) / 
	(SolnTime[iFl] - SolnTime[iPr]);
      Wtp = 1.0 - Wtf;

      /* Must be after the final Soln */
    } else {
      Wtp = 0.0;
      Wtf = 1.0;      
    }

    /* Apply Soln data to Cal - deal with bad solutions */
    if ((outCalRow->dAz!=fblank) && (SolndAz[iPr]!=fblank) && (SolndAz[iFl]!=fblank))
      outCalRow->dAz += (Wtp * SolndAz[iPr]  + Wtf * SolndAz[iFl]);
    else
      outCalRow->dAz = fblank;

    if ((outCalRow->dEl!=fblank) && (SolndEl[iPr]!=fblank) && (SolndEl[iFl]!=fblank))
      outCalRow->dEl += (Wtp * SolndEl[iPr] + Wtf * SolndEl[iFl]);
    else
      outCalRow->dEl = fblank;

     /* Additive term */
    for (i=0; i<ndetect; i++) {
      if ((outCalRow->add[i]!=fblank) && (SolnAdd[iPr][i]!=fblank) && (SolnAdd[iFl][i]!=fblank)) {
	if (fabs(outCalRow->mult[i]) > 1.0e-10) minMult= outCalRow->mult[i];
	else minMult = 1.0e-10;
	outCalRow->add[i]  +=  (Wtp * SolnAdd[iPr][i]+ Wtf * SolnAdd[iFl][i]) / minMult;
      }
      else
	outCalRow->add[i] = fblank;
      /* Debug Mack Hack Hack
	 if (outCalRow->add[i] != fblank) {
	 outCalRow->add[i] = MIN (outCalRow->add[i],  0.2);
	 outCalRow->add[i] = MAX (outCalRow->add[i], -0.2);
	 } */
    }

    /* Cal value */
    for (i=0; i<ndetect; i++) {
      if ((outCalRow->cal[i]!=fblank) && (SolnCal[iPr][i]!=fblank) && (SolnCal[iFl][i]!=fblank)) {
	if (fabs(outCalRow->mult[i]) > 1.0e-10) minMult= outCalRow->mult[i];
	else minMult = 1.0e-10;
	outCalRow->cal[i]  +=  (Wtp * SolnCal[iPr][i]+ Wtf * SolnCal[iFl][i]) / minMult;
      }
      else {
	outCalRow->cal[i] = fblank;
	outCalRow->add[i] = fblank;
      }
    }
    
    /* Polynomial term */
    for (i=0; i<npoly; i++) {
      if ((outCalRow->poly[i]!=fblank) && (SolnPoly[iPr][i]!=fblank) && (SolnPoly[iFl][i]!=fblank)) {
	outCalRow->poly[i] += (Wtp * SolnPoly[iPr][i] + Wtf * SolnPoly[iFl][i]);
      } else {
	outCalRow->poly[i] = fblank;
      }
    }

    /* Multiplicative term */
     for (i=0; i<ndetect; i++) {
      if ((outCalRow->mult[i]!=fblank) && (SolnMult[iPr][i]!=fblank) && (SolnMult[iFl][i]!=fblank))
	outCalRow->mult[i] *= (Wtp * SolnMult[iPr][i]+ Wtf * SolnMult[iFl][i]);
      else {
	outCalRow->mult[i] = fblank;
 	outCalRow->add[i] = fblank;
      }
   }

    /* Weight */
     for (i=0; i<ndetect; i++) {
      if ((outCalRow->mult[i]!=fblank) && (SolnWt[iPr][i]!=fblank) && (SolnWt[iFl][i]!=fblank))
	outCalRow->wt[i] *= (Wtp * SolnWt[iPr][i]+ Wtf * SolnWt[iFl][i]);
      else {
	outCalRow->wt[i] = fblank;
 	outCalRow->add[i] = fblank;
      }
   }

     /* debug
	if ((outCalRow->Time>2.0210e-4)&&(outCalRow->Time<2.0250e-4)) {
	fprintf (stderr, "S2C: time %e prior %e %e follow %e %e\n",
	outCalRow->Time, SolnTime[iPr], Wtp, SolnTime[iFl], Wtf);
	} */
     
   outCalRow->status = 1;

    /* write row */
    retCode = ObitTableOTFCalWriteRow (outCal, iRow, outCalRow, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outCal);
  } /* end loop calibrating Cal table */
  
 
 /* Close cal tables */
  if ((ObitTableOTFSolnClose (inSoln, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input OTFSoln Table file");
    return outCal;
  }
  
  if ((ObitTableOTFCalClose (inCal, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input OTFCal Table file");
    return outCal;
  }
  
  if ((ObitTableOTFCalClose (outCal, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output OTFCal Table file");
    return outCal;
  }
  
  /* Cleanup */
  inSolnRow = ObitTableOTFSolnRowUnref (inSolnRow);
  inCalRow  = ObitTableOTFCalRowUnref (inCalRow);
  outCalRow = ObitTableOTFCalRowUnref (outCalRow);
  if (SolnCal[0])  g_free(SolnCal[0]);
  if (SolnAdd[0])  g_free(SolnAdd[0]);
  if (SolnMult[0]) g_free(SolnMult[0]);
  if (SolnWt[0])   g_free(SolnWt[0]);
  if (SolnPoly[0]) g_free(SolnPoly[0]);
  if (SolnCal[1])  g_free(SolnCal[1]);
  if (SolnAdd[1])  g_free(SolnAdd[1]);
  if (SolnMult[1]) g_free(SolnMult[1]);
  if (SolnWt[1])   g_free(SolnWt[1]);
  if (SolnPoly[1]) g_free(SolnPoly[1]);
 

  return outCal;
} /* end ObitOTFSoln2Cal */

/*----------------------Private functions---------------------------*/
/**
 * Copy an open Soln table to an open Cal table
 * \param inSoln   Input Soln table
 * \param outCal   Output Cal table.
 * \param err      Error stack, returns if not empty.
 */
void ObitOTFSolnCopyCal (ObitTableOTFSoln *inSoln, ObitTableOTFCal *outCal, 
			 ObitErr *err)
{
  olong i, irow, limit, npoly, ndetect;
  ObitIOCode retCode;
  ObitTableOTFSolnRow *SolnRow=NULL;
  ObitTableOTFCalRow  *outCalRow=NULL;
  gchar *routine = "ObitOTFSolnCopyCal";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitTableOTFSolnIsA(inSoln));
  g_assert (ObitTableOTFCalIsA(outCal));

  /* Create Soln Row */
  SolnRow = newObitTableOTFSolnRow (inSoln);

  /* Attach Soln row to output buffer */
  ObitTableOTFSolnSetRow (inSoln, SolnRow, err);
  if (err->error) Obit_traceback_msg (err, routine, inSoln->name);

  /* Create output Row */
  outCalRow = newObitTableOTFCalRow (outCal);

  /* Attach row to output buffer */
  ObitTableOTFCalSetRow (outCal, outCalRow, err);
  if (err->error) Obit_traceback_msg (err, routine, inSoln->name);

  /* how many of things? */
  npoly   = inSoln->numPoly;
  ndetect = inSoln->numDet;

  /* loop  converting Soln to Cal table */
  retCode = OBIT_IO_OK;
  limit = inSoln->myDesc->nrow; /* how many rows to read */
  for (irow=0; irow<=limit; irow++) {

    /* read row */
    retCode = ObitTableOTFSolnReadRow (inSoln, irow, SolnRow, err);
    if (err->error) Obit_traceback_msg (err, routine, inSoln->name);
    if (retCode==OBIT_IO_EOF) break; /* done? */

    /* Copy Soln data to Cal */
    outCalRow->Time = SolnRow->Time;
    outCalRow->TimeI = SolnRow->TimeI;
    outCalRow->Target = SolnRow->Target;
    outCalRow->dAz = SolnRow->dAz;
    outCalRow->dEl = SolnRow->dEl;
    for (i=0; i<ndetect; i++) outCalRow->cal[i]  = SolnRow->cal[i];
    for (i=0; i<ndetect; i++) outCalRow->add[i]  = SolnRow->add[i];
    for (i=0; i<ndetect; i++) outCalRow->mult[i] = SolnRow->mult[i];
    for (i=0; i<ndetect; i++) outCalRow->wt[i]   = SolnRow->wt[i];
    for (i=0; i<npoly; i++)   outCalRow->poly[i] = SolnRow->poly[i];
    outCalRow->status = 1;

    /* write row */
    retCode = ObitTableOTFCalWriteRow (outCal, irow, outCalRow, err);
    if (err->error) Obit_traceback_msg (err, routine, inSoln->name);
  } /* end loop converting Soln to Cal table */
  
 /* Close cal tables */
  if ((ObitTableOTFSolnClose (inSoln, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing input OTFSoln Table file");
    return;
  }
  
  if ((ObitTableOTFCalClose (outCal, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output OTFCal Table file");
    return;
  }
  
  /* Cleanup */
  SolnRow = ObitTableOTFSolnRowUnref (SolnRow);
  outCalRow = ObitTableOTFCalRowUnref (outCalRow);

} /* end ObitOTFSolnCopyCal */
