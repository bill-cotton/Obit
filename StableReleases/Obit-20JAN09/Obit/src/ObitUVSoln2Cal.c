/* $Id$   */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitUVSoln2Cal.h"
#include "ObitTableUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVSoln2Cal.c
 * ObitUVSoln2Cal utility calibration function definitions.
 */

/*---------------Private function prototypes----------------*/
/* Copy a Soln table to a Cal table */
void ObitUVSolnCopyCal (ObitTableSN *inSoln, ObitTableCL *outCal,  ObitErr *err);

/*----------------------Public functions---------------------------*/
/**
 * Apply a Soln (SN) table to a Cal (CL) table and write a new Cal table
 * Calibration parameters are on the inUV info member.
 * If an input Cal table is given then apply Solutions in this routine,
 * if no input Cal table, then copy the Soln table to a new Cal table
 * in ObitUVSolnCopyCal.
 * Input SN table will have its phases referenced to refAnt.
 * \param inUV    Input UV data. 
 * Control parameters on inUV:
 * \li  "Sources" OBIT_string (?,?,1) Source names selected unless any starts with
 *                 a '-' in which case all are deselected (with '-' stripped).
 * \li  "souCode"  OBIT_string (4,1,1) Source Cal code desired, '    ' => any code selected
 *                                   '*   ' => any non blank code (calibrators only)
 *                                   '-CAL' => blank codes only (no calibrators)
 * \li  "Qual"     Obit_int (1,1,1)  Source qualifier, -1 [default] = any
 * \li  "calSour"  OBIT_string (?,?,1) Calibrator names selected unless any 
 *                 starts with a '-' in which cse all are deselected 
 *                 (with '-' stripped).
 * \li "calCode"   OBIT_strine (4,1,1) Calibrator code
 * \li "solnVer"   OBIT_int (1,1,1) Input Solution (SN) table version 
 * \li "calIn"     OBIT_int (1,1,1) Input Cal (CL) table version 
 *                 iff <0 then no input Cal table, copy Soln records to output.
 * \li "calOut"    OBIT_int (1,1,1) Output Calibration table version
 * \li  "Antennas" OBIT_int (?,1,1) a list of selected antenna numbers, if any is negative
 *                 then the absolute values are used and the specified antennas are deselected.
 * \li "subA"      OBIT_int (1,1,1) Selected subarray (default 1)
 * \li "refAnt"    OBIT_int (1,1,1) Ref ant to use. (default 1)
 * \li "allPass", OBIT_bool (1,1,1) If true copy unmodified entries as well.
 *                else only data calibrated.  Default = FALSE.
 * \param outUV   UV with which the output  UVCal is to be associated
 * \param err     Error stack, returns if not empty.
 * \return Pointer to the newly created CL object which is associated with outUV.
 *         Should be Unreffed when done.
 */
ObitTableCL* ObitUVSoln2Cal (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitTableCL    *outCal=NULL, *inCal=NULL;
  ObitTableSN    *inSoln=NULL;
  ObitTableSNRow *SolnRow=NULL;
  ObitTableCLRow *CalRow=NULL;
  ObitUVSoln *mySoln = NULL;
  ObitUVSel *sel;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gboolean want, allPass, warnMaxInter;
  olong iRow, oRow, ver, highVer, limitC;
  olong  solnVer, calIn, calOut, refAnt, subA, iif,itemp;
  oint numPol, numIF, numTerm;
  ObitIOCode retCode;
  ofloat tr, ti, gr, gi;
  ofloat fblank = ObitMagicF();
  gchar *tname;
  gchar        *solnParms[] = {  /* Control Parameters to copy to ObitUVSoln */
    "interMode", "interParm", "interNPoly", "maxInter", "solnVer", 
    "doBlank", "allPass", "smoFunc", 
    "Sources", "souCode", "Qual", "calSour", "calCode",
    "FreqID", "timeRange", "subA", "Antennas", "refAnt", 
    "solnVer", "calIn", "calOut", 
    NULL
  };
  gchar *routine = "ObitUVSoln2Cal";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outCal;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));

  /* Get table numbers and other parameters */
  solnVer = 0;
  ObitInfoListGetTest(inUV->info, "solnVer", &type, dim, &solnVer);
  calIn = 0;
  ObitInfoListGetTest(inUV->info, "calIn", &type, dim, &calIn);
  calOut = 0;
  ObitInfoListGetTest(inUV->info, "calOut", &type, dim, &calOut);
  /* Refuse to write to CL 1 */
  Obit_retval_if_fail((calOut!=1), err, outCal,
		      "%s: You may not modify CL table 1", routine);
  subA = 1;
  ObitInfoListGetTest(inUV->info, "subA", &type, dim, &subA);
  refAnt = 1;
  ObitInfoListGetTest(inUV->info, "refAnt", &type, dim, &refAnt);
  itemp = 0;
  ObitInfoListGetTest(inUV->info, "solnVer", &type, dim, &itemp);
  solnVer = itemp;
  allPass = FALSE;
  ObitInfoListGetTest(inUV->info, "allPass", &type, dim, &allPass);

   /* Full instantiation of uv object */
  ObitUVFullInstantiate (inUV, TRUE, err);
  if (err->error) goto cleanup;

  /* Which SN table? */
  highVer = ObitTableListGetHigh (inUV->tableList, "AIPS SN");
  if (solnVer==0) solnVer = highVer;
  ver = solnVer;
  inSoln = newObitTableSNValue (inUV->name, (ObitData*)inUV, &ver, 
				OBIT_IO_ReadOnly, 0, 0, err);
  if (err->error) goto cleanup;

  /* Check SN table */
  if (!ObitTableSNIsA(inSoln)) {
    Obit_log_error(err, OBIT_Error, "%s: SN table %d not found in %s",
		   routine, solnVer, inUV->name);
    goto cleanup;
  }

  /* Sort input SN if needed */
  ObitTableUtilSort ((ObitTable*)inSoln, "TIME    ", FALSE, err);
  if (err->error) goto cleanup;

  /* Rereference SN table */
  ObitUVSolnRefAnt (inSoln, subA, &refAnt, err);
  if (err->error) goto cleanup;

  /* if input Cal table exists then open */
  if (calIn>=0) {
    /* Instantiate input calibration table object if one exists */
    tname = g_strconcat ("Input Calibration for: ",inUV->name, NULL);
    /* Which CL table? */
    ver = calIn;
    highVer = ObitTableListGetHigh (inUV->tableList, "AIPS CL");
    if (ver==0) ver = highVer;
    inCal = newObitTableCLValue(tname, (ObitData*)outUV, &ver, 
				OBIT_IO_ReadWrite,  0, 0, 0, err);
    g_free (tname);
    if (err->error) goto cleanup;
    
    /* Check CL table */
    if (!ObitTableCLIsA(inCal)) {
      Obit_log_error(err, OBIT_Error, "%s: CL table %d not found in %s",
		     routine, calIn, inUV->name);
      goto cleanup;
    }

    /* Sort input CL if needed */
    ObitTableUtilSort ((ObitTable*)inCal, "TIME    ", FALSE, err);
    if (err->error) goto cleanup;

    /* Open input table */
    ObitTableCLOpen (inCal, OBIT_IO_ReadWrite, err);
    if (err->error) goto cleanup;
    
    /* Better have some data */
    if (inCal->myDesc->nrow<=0) {
      Obit_log_error(err, OBIT_Error, "Input CL table has no data");
      goto cleanup;
    }

  } /* end open input Cal table */

  /* Instantiate/Create output */
  tname = g_strconcat ("Output Calibration for: ",outUV->name, NULL);
  if (calIn>=0) { /* input Cal table? */
    numPol  = inCal->numPol;
    numIF   = inCal->numIF;
    numTerm = inCal->numTerm;
  } else { /* Nope - only Soln */
    numPol  = inSoln->numPol;
    numIF   = inSoln->numIF;
    numTerm = 1;
  }
  /* Which CL table? */
  ver = calOut;
  highVer = ObitTableListGetHigh (outUV->tableList, "AIPS CL");
  if (ver==0) ver = highVer+1;
  outCal = newObitTableCLValue(tname, (ObitData*)outUV, &ver, OBIT_IO_ReadWrite,  
			       numPol, numIF, numTerm, err);
  g_free (tname);
  if (err->error) goto cleanup;

  /* if no input Cal table exists then copy input Soln table to output Cal table */
  if (calIn<0) {
    ObitUVSolnCopyCal (inSoln, outCal, err);
    goto cleanup;
 }

  /* Inform user */
  Obit_log_error(err, OBIT_InfoErr, 
			 "Apply SN table %d to CL table %d and write CL table %d", 
			 solnVer, calIn, calOut);

  /* Create solution interpolator */
  mySoln = ObitUVSolnCreate ("Soln interpolator", inUV);

  /* Copy control parameters */
  ObitInfoListCopyList (inUV->info, mySoln->info, solnParms);
  
  /* Startup */
  ObitUVSolnStartUp (mySoln, err);
  if (err->error) goto cleanup;

  /* Startup may have smoothed inSoln to another table */
  if (inSoln!=mySoln->SNTable) {
    inSoln = ObitTableSNUnref(inSoln);
    inSoln = ObitTableSNRef(mySoln->SNTable);
  }

  /* Create work Soln Row */
  SolnRow = newObitTableSNRow (inSoln);

  /* Attach row to buffer to allow for storage - does this cause trouble? */
  ObitTableSNSetRow (inSoln, SolnRow, err);

  /* Deselect any prior entries to be replaced 
     this needs to be after ObitUVSolnStartUp call so selector is set */
  sel = inUV->mySel;
  ObitUVSolnDeselCL (outCal,  sel->SubA, sel->FreqID, 
		     sel->numberAntList, sel->ants, 
		     sel->numberSourcesList, sel->sources, 
		     sel->timeRange, err);
  if (err->error) goto cleanup;

  /* Open output table */
  ObitTableCLOpen (outCal, OBIT_IO_ReadWrite, err) ;
  if (err->error) goto cleanup;

  /* Create CL table Row */
  CalRow = newObitTableCLRow (outCal);

  /* Attach row to output buffer */
  ObitTableCLSetRow (outCal, CalRow, err);
  if (err->error) goto cleanup;

  /* Are there previous entries in outCal - 
     if so output probably not in time order when we get through here */
  if (outCal->myDesc->nrow>0) outCal->myDesc->sort[0] = -1;
  /* else it will be time ordered */
  else {
    outCal->myDesc->sort[0] = outCal->TimeCol + 1;
    outCal->numAnt = inCal->numAnt;
  }
  
  warnMaxInter = FALSE;
  retCode = OBIT_IO_OK;
  limitC = inCal->myDesc->nrow;  /* how many Cal rows to read */

  /* Loop interpolating Soln table and applying to input Cal table */
  for (iRow=1; iRow<=limitC; iRow++) {

    /* read input row */
    retCode = ObitTableCLReadRow (inCal, iRow, CalRow, err);
    if (err->error) goto cleanup;
    if (retCode==OBIT_IO_EOF) break; /* done? */
    if (CalRow->status < 0) continue; /* entry flagged? */

    /* Set selection */
    SolnRow->Time   = CalRow->Time;
    SolnRow->TimeI  = CalRow->TimeI;
    SolnRow->SourID = CalRow->SourID;
    SolnRow->antNo  = CalRow->antNo;
    SolnRow->SubA   = CalRow->SubA;
    SolnRow->FreqID = CalRow->FreqID;

    /* Interpolate calibration to SolnRow */
    want = ObitUVSolnGetSN (mySoln, SolnRow, err);
    if (err->error) goto cleanup;
    /* Not to modify but want to copy original? */
    if (!want && (allPass)) goto writeit;
    if (!want) continue;  /* This one selected? */

    /* Check if cal entry out of maxInter range */
    warnMaxInter = warnMaxInter ||
      (mySoln->PriorCalTime-mySoln->maxInter > mySoln->CalTime) ||
      (mySoln->FollowCalTime+mySoln->maxInter < mySoln->CalTime) ||
      ((mySoln->PriorCalTime+mySoln->maxInter < mySoln->CalTime) &&
       (mySoln->FollowCalTime-mySoln->maxInter > mySoln->CalTime));

    /* Apply calibration  */
    CalRow->MBDelay1 += SolnRow->MBDelay1;
    for (iif=0; iif<numIF; iif++) {
      if ((SolnRow->Weight1[iif]>0.0) && (CalRow->Real1[iif]!=fblank)) {
	tr = CalRow->Real1[iif];
	ti = CalRow->Imag1[iif];
	gr = SolnRow->Real1[iif];
	gi = SolnRow->Imag1[iif];
	CalRow->Real1[iif]    = tr*gr - ti*gi;
	CalRow->Imag1[iif]    = tr*gi + ti*gr;
	CalRow->Rate1[iif]   += SolnRow->Rate1[iif];
	CalRow->Delay1[iif]  += SolnRow->Delay1[iif];
	CalRow->Weight1[iif] += SolnRow->Weight1[iif];
	CalRow->RefAnt1[iif]  = SolnRow->RefAnt1[iif];
      } else {
	CalRow->Real1[iif]  = fblank;
	CalRow->Imag1[iif]  = fblank;
	CalRow->Rate1[iif]  = fblank;
	CalRow->Delay1[iif] = fblank;
	CalRow->Weight1[iif]= 0.0;
      }
    }
    if (numPol>1) {
      CalRow->MBDelay2 += SolnRow->MBDelay2;
      for (iif=0; iif<numIF; iif++) {
	if ((SolnRow->Weight2[iif]>0.0) && (CalRow->Real2[iif]!=fblank)) {
	  tr = CalRow->Real2[iif];
	  ti = CalRow->Imag2[iif];
	  gr = SolnRow->Real2[iif];
	  gi = SolnRow->Imag2[iif];
	  CalRow->Real2[iif]    = tr*gr - ti*gi;
	  CalRow->Imag2[iif]    = tr*gi + ti*gr;
	  CalRow->Rate2[iif]   += SolnRow->Rate2[iif];
	  CalRow->Delay2[iif]  += SolnRow->Delay2[iif];
	  CalRow->Weight2[iif] += SolnRow->Weight2[iif];
	  CalRow->RefAnt2[iif]  = SolnRow->RefAnt2[iif];
	} else {
	  CalRow->Real2[iif]  = fblank;
	  CalRow->Imag2[iif]  = fblank;
	  CalRow->Rate2[iif]  = fblank;
	  CalRow->Delay2[iif] = fblank;
	  CalRow->Weight2[iif]= 0.0;
	}
      }
    }

    /* write row */
  writeit:
    oRow = -1;
    retCode = ObitTableCLWriteRow (outCal, oRow, CalRow, err);
    if (err->error) goto cleanup;
  } /* end loop calibrating Cal table */
  
  if (outCal) ObitTableCLClose (outCal, err);

  /* Sort output CL if needed */
  ObitTableUtilSort ((ObitTable*)outCal, "TIME    ", FALSE, err);
  if (err->error) goto cleanup;

  /* Check if time lost to maxInter */
  if (warnMaxInter) {
    Obit_log_error(err, OBIT_InfoWarn, 
		   "Output CL entries flagged due to maxInter");
  }
  
  /* Cleanup */
 cleanup:
  /* Close cal tables */
  if (inCal)  ObitTableCLClose (inCal, err);
  if (outCal) ObitTableCLClose (outCal, err);
  if (mySoln) ObitUVSolnShutDown (mySoln, err); 
  mySoln    = ObitUVSolnUnref (mySoln);
  inSoln    = ObitTableSNUnref(inSoln);
  inCal     = ObitTableCLUnref(inCal);
  SolnRow   = ObitTableSNRowUnref (SolnRow);
  CalRow = ObitTableCLRowUnref (CalRow);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outCal);
 
  return outCal;
} /* end ObitUVSoln2Cal */

/*----------------------Private functions---------------------------*/
/**
 * Copy a Soln (SN) table to a Cal (CL) table
 * \param inSoln  Input Solution (SN) table
 * \param outCal   Output Cal table.
 * \param err      Error stack, returns if not empty.
 */
void ObitUVSolnCopyCal (ObitTableSN *inSoln, ObitTableCL *outCal, ObitErr *err)
{
  olong i, irow, orow, numPol, numIF, limit;
  ObitIOCode retCode;
  ObitTableSNRow  *SolnRow=NULL;
  ObitTableCLRow  *CalRow=NULL;
  gchar *routine = "ObitUVSolnCopyCal";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitTableSNIsA(inSoln));
  g_assert (ObitTableCLIsA(outCal));

  /* Open SN table */
  ObitTableSNOpen (inSoln, OBIT_IO_ReadOnly, err);
  if (err->error) goto cleanup;

  /* Create Soln Row */
  SolnRow = newObitTableSNRow (inSoln);

  /* Open output table */
  ObitTableCLOpen (outCal, OBIT_IO_ReadWrite, err) ;
  if (err->error) goto cleanup;

  /* Create output CL table Row */
  CalRow = newObitTableCLRow (outCal);

  /* Attach row to output buffer */
  ObitTableCLSetRow (outCal, CalRow, err);
  if (err->error) goto cleanup;

  /* how many of things? */
  numIF  = inSoln->numIF;
  numPol = inSoln->numPol;

  /* loop  converting Soln to Cal table */
  retCode = OBIT_IO_OK;
  limit = inSoln->myDesc->nrow; /* how many rows to read */
  for (irow=0; irow<=limit; irow++) {

    /* read row */
    retCode = ObitTableSNReadRow (inSoln, irow, SolnRow, err);
    if (err->error) Obit_traceback_msg (err, routine, inSoln->name);
    if (retCode==OBIT_IO_EOF) break; /* done? */

    /* Copy Soln data to Cal */
    CalRow->Time   = SolnRow->Time;
    CalRow->TimeI  = SolnRow->TimeI;
    CalRow->SourID = SolnRow->SourID;
    CalRow->antNo  = SolnRow->antNo;
    CalRow->SubA   = SolnRow->SubA;
    CalRow->FreqID = SolnRow->FreqID;
    CalRow->IFR    = SolnRow->IFR;
    CalRow->atmos  = 0.0;
    CalRow->Datmos = 0.0;
    for (i=0; i<outCal->myDesc->repeat[outCal->GeoDelayCol]; i++) 
      CalRow->GeoDelay[i] = 0.0;
    for (i=0; i<numIF; i++) CalRow->DopplerOff[i] = 0.0;
    CalRow->MBDelay1  = SolnRow->MBDelay1;
    CalRow->clock1    = 0.0;
    CalRow->Dclock1   = 0.0;
    CalRow->dispers1  = 0.0;
    CalRow->Ddispers1 = 0.0;
    for (i=0; i<numIF; i++) CalRow->Real1[i]  = SolnRow->Real1[i];
    for (i=0; i<numIF; i++) CalRow->Imag1[i]  = SolnRow->Imag1[i];
    for (i=0; i<numIF; i++) CalRow->Rate1[i]  = SolnRow->Rate1[i];
    for (i=0; i<numIF; i++) CalRow->Delay1[i] = SolnRow->Delay1[i];
    for (i=0; i<numIF; i++) CalRow->Weight1[i]= SolnRow->Weight1[i];
    for (i=0; i<numIF; i++) CalRow->RefAnt1[i]= SolnRow->RefAnt1[i];
    if (numPol>1) {
      CalRow->MBDelay2  = SolnRow->MBDelay2;
      CalRow->clock2    = 0.0;
      CalRow->Dclock2   = 0.0;
      CalRow->dispers2  = 0.0;
      CalRow->Ddispers2 = 0.0;
      for (i=0; i<numIF; i++) CalRow->Real2[i]  = SolnRow->Real2[i];
      for (i=0; i<numIF; i++) CalRow->Imag2[i]  = SolnRow->Imag2[i];
      for (i=0; i<numIF; i++) CalRow->Rate2[i]  = SolnRow->Rate2[i];
      for (i=0; i<numIF; i++) CalRow->Delay2[i] = SolnRow->Delay2[i];
      for (i=0; i<numIF; i++) CalRow->Weight2[i]= SolnRow->Weight2[i];
      for (i=0; i<numIF; i++) CalRow->RefAnt2[i]= SolnRow->RefAnt2[i];
    }

    CalRow->status = 1;

    /* write row */
    orow = -1;
    retCode = ObitTableCLWriteRow (outCal, orow, CalRow, err);
    if (err->error) goto cleanup;
  } /* end loop converting Soln to Cal table */
  
  /* Cleanup */
 cleanup:
  ObitTableSNClose (inSoln, err);
  ObitTableCLClose (outCal, err);
  SolnRow = ObitTableSNRowUnref (SolnRow);
  CalRow  = ObitTableCLRowUnref (CalRow);
  if (err->error) Obit_traceback_msg (err, routine, inSoln->name);

} /* end ObitUVSolnCopyCal */
