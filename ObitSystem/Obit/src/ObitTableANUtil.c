/* $Id$ */
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

#include <math.h>
#include "ObitTableANUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableANUtil.c
 * ObitTableAN class utility function definitions.
 */

/*----------------------Public functions---------------------------*/


/**
 * Reads and returns information from an "AIPS AN" table.
 * Undesired output parameters are indicated by NULL.
 * \param in       Table to obtain data from
 * \param numAnt   Highest antenna number if non NULL.
 * \param refFreq  Reference frequency if non NULL
 * \param refDate  Reference data if non NULL, must be at least 11 characters allocated.
 * \param *err     ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableANGetInfo (ObitTableAN *in, oint *numAnt, odouble *refFreq, 
			       gchar* refDate, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableANRow *row;
  oint maxAnt = 0;
  gchar *routine = "ObitTableANGetInfo";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitTableANIsA(in));

  /* Open table */
  retCode = ObitTableANOpen (in, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine,in->name, retCode);

  /* Get items from the header */
  if (refFreq!=NULL) *refFreq = in->Freq;
  if (refDate!=NULL) strncpy (refDate, in->RefDate, 10);

  /* Only need to read table if numAnt != NULL */
  if (numAnt != NULL) {

    /* Create table row */
    row = newObitTableANRow (in);

    /* loop over table */
    while (retCode==OBIT_IO_OK) {
      retCode = ObitTableANReadRow (in, -1, row, err);
      if (retCode == OBIT_IO_EOF) break;
      
      /* Get maximum antenna number */
      maxAnt = MAX (maxAnt, row->noSta);
    }
    
    /* check for errors */
    if ((retCode > OBIT_IO_EOF) || (err->error))
      Obit_traceback_val (err, routine,in->name, retCode);
    
    /* Release table row */
    row = ObitTableANRowUnref (row);

    /* return maximum antenna number */
    *numAnt = maxAnt;

  } /* end of reading table */

    /* Close table */
  retCode = ObitTableANClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine,in->name, retCode);
 

  return retCode;
} /* end ObitTableANGetInfo */

/**
 * Convert the contents of a ObitTableAN into an ObitAntennaList.
 * \param in    Table to obtain data from
 * \param *err  ObitErr error stack.
 * \return requested ObitAntennaList
 */
ObitAntennaList* ObitTableANGetList (ObitTableAN *in, ObitErr *err) {
  ObitAntennaList *out=NULL;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableANRow *row;
  olong irow;
  olong maxANid, i, iant;
  ObitInfoType type;
  gboolean doVLA, doVLBI, doATCA;
  odouble x, y, z, ArrLong, rho, dtemp;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar tempName[101]; /* should always be big enough */
  gchar *routine = "ObitTableANGetList";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitTableANIsA(in));

  /* Open table */
  retCode = ObitTableANOpen (in, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, out);
  
  /* Create table row */
  row = newObitTableANRow (in);

  /* loop over table looking for highest number */
  maxANid = -1;
  irow = 0;
  while (retCode==OBIT_IO_OK) {
    irow++;
    retCode = ObitTableANReadRow (in, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, out);
    
    maxANid = MAX (maxANid, row->noSta);

  } /* end loop over file */
  
  /* check for errors */
  if ((retCode > OBIT_IO_EOF) || (err->error))
    Obit_traceback_val (err, routine, in->name, out);

  /* Create output */
  g_snprintf (tempName, 100, "Antenna List for %s",in->name);
  out = ObitAntennaListCreate (tempName, maxANid, in->numPCal);
  
  /* Get table header information */
  /* Have to see of polarization type is in the header ('POLTYPE') */
  g_snprintf (tempName, 7, "       ");
  dim[0] = 7; dim[1] = dim[2] = dim[3] = 1;;
  ObitInfoListGetTest(in->myDesc->info, "POLTYPE", &type, dim, &tempName);
  tempName[dim[0]] = 0;

  out->polType     = ObitAntennaListGetPolType(tempName);
  ObitUVDescDate2JD (in->RefDate, &out->JD);
  out->ArrayXYZ[0] = in->ArrayX;
  out->ArrayXYZ[1] = in->ArrayY;
  out->ArrayXYZ[2] = in->ArrayZ;
  out->GSTIAT0     = in->GSTiat0*DG2RAD;   /* Convert to radians */
  out->RotRate     = in->DegDay*DG2RAD;    /* Convert to radians/day */
  out->PolarXY[0]  = in->PolarX;
  out->PolarXY[1]  = in->PolarY;
  out->ut1Utc      = in->ut1Utc*DG2RAD*15./3600.0;  /* Convert to radians */
  out->dataUtc     = in->dataUtc*DG2RAD*15./3600.0; /* Convert to radians */
  if (!strncmp (in->TimeSys, "IAT", 3))
    out->dataIat = 0.0;
  else
    out->dataIat  = in->iatUtc*DG2RAD*15./3600.0; /* Convert to radians */
  out->numPoln     = in->numPCal;
  for (i=0; i<12; i++) out->TimSys[i]  = in->TimeSys[i];
  for (i=0; i<12; i++) out->ArrName[i] = in->ArrName[i];
  out->FreqID      = in->FreqID;
  out->numPCal = in->numPCal;
  out->numIF   = in->numPCal/2; /* look out */

  /* Polarization reference antenna */
  out->polRefAnt = 1;
  ObitInfoListGetTest(in->myDesc->info, "P_REFANT", &type, dim, &out->polRefAnt);

  /* Polarization R-L Phases */
  out->RLPhaseDiff = g_realloc (out->RLPhaseDiff, out->numIF*sizeof(ofloat));
  for (i=0; i<out->numIF; i++) {
    g_snprintf (tempName, 9, "P_DIFF%2.2d",i+1);
    dtemp = 0.0;
    ObitInfoListGetTest(in->myDesc->info, tempName, &type, dim, &dtemp);
    out->RLPhaseDiff[i] = dtemp;
  } 

  /* Some array dependent information */
  /* Is this the VLA? */
  doVLA  = !strncmp(in->ArrName, "VLA     ", 8);
  out->isVLA = doVLA;

  /* Is this the ATCA? It uses earth center but without Y flip like VLBI */
  doATCA = !strncmp(in->ArrName, "ATCA    ", 8);

  /* Otherwise VLBI Uses earth center, but with Y with sign flip */
  doVLBI = (!doATCA ) && 
    (fabs(in->ArrayX)<1000.0) && (fabs(in->ArrayY)<1000.0) && (fabs(in->ArrayZ)<1000.0);

  /* loop over table saving information */
  retCode = OBIT_IO_OK;
  irow = 0;
  while (retCode==OBIT_IO_OK) {
    irow++;
    retCode = ObitTableANReadRow (in, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, out);

    iant = row->noSta-1;
    /* Check antenna number */
    Obit_retval_if_fail((iant>=0), err, out,
			"%s: Corrupt antenna table on %s", 
			routine, in->name);

    out->ANlist[iant]->AntID     = row->noSta;
    out->ANlist[iant]->AntMount  = row->mntSta;
    out->ANlist[iant]->numPCal   = out->numPoln;
    out->ANlist[iant]->FeedAPA   = row->PolAngA;
    out->ANlist[iant]->FeedBPA   = row->PolAngB;
    out->ANlist[iant]->FeedAType  = row->polTypeA[0];
    out->ANlist[iant]->FeedBType  = row->polTypeB[0];
    for (i=0; i<8; i++) out->ANlist[iant]->AntName[i]  = row->AntName[i];
    for (i=0; i<3; i++) out->ANlist[iant]->AntXYZ[i]  = row->StaXYZ[i];
    for (i=0; i<out->numPoln; i++) out->ANlist[iant]->FeedAPCal[i]  = row->PolCalA[i];
    for (i=0; i<out->numPoln; i++) out->ANlist[iant]->FeedBPCal[i]  = row->PolCalB[i];

    /* Lat, (E) long in radians
       X => (lat 0, long 0)
       Y => (lat 0, long 90E)
       Z => (lat 90) */

    x = in->ArrayX;
    y = in->ArrayY;
    z = in->ArrayZ;
    
   
    /* Adjust coordinates  by telescope */
    if (doVLA) {
      ArrLong = 1.878283678;
      x += row->StaXYZ[0]*cos (ArrLong) + row->StaXYZ[1]*sin (ArrLong);
      y += row->StaXYZ[1]*cos (ArrLong) - row->StaXYZ[0]*sin (ArrLong);
      z += row->StaXYZ[2];
    } else if (doVLBI) {
      x = x + row->StaXYZ[0];
      y = -(y + row->StaXYZ[1]);  /* Flip handedness of VLBI data */
      z = z + row->StaXYZ[2];
    } else if (doATCA) {
      x = x + row->StaXYZ[0];
      y = y + row->StaXYZ[1];
      z = z + row->StaXYZ[2];
    }
    
    /* Get lat/long/radius */
    rho = sqrt(x*x + y*y + z*z);
    if (x!=0.0) out->ANlist[iant]->AntLong = atan2 (y, x);
    else out->ANlist[iant]->AntLong = 0.0;
    if(rho!=0.0) out->ANlist[iant]->AntLat  = asin(z/rho);
    else out->ANlist[iant]->AntLat  = 0.0;
    out->ANlist[iant]->AntRad  = rho;

    /* Make sure first antenna filled in (used for polarization cal ) */
    if ((iant>0) && (out->ANlist[0]<=0)) {
      out->ANlist[0]->AntLong = out->ANlist[iant]->AntLong;
      out->ANlist[0]->AntLat  = out->ANlist[iant]->AntLat;
      out->ANlist[0]->AntRad  = out->ANlist[iant]->AntRad;      
    }

  } /* end second loop over table */

 /* Release table row */
  row = ObitTableANRowUnref (row);

  /* Close table */
  retCode = ObitTableANClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, out);

  return out;
} /* end ObitTableANGetList */

/**
 * Copies AN tables from inUV to outUV with selection in inUV
 * If poln calibration is selected on inUV, polarization 
 * calibration info is removed.
 * \param inUV     Input UV to copy from
 * \param outUV    Output UV to copy to
 * \param *err     ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableANSelect (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableAN    *inTab=NULL, *outTab=NULL;
  ObitTableANRow *inRow=NULL, *outRow=NULL;
  ObitErr *terr=NULL;
  ObitInfoType type;
  olong iif, oif, i, polRefAnt;
  olong highANver, iANver, inANRow, outANRow;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  oint numOrb, numIF, numPCal;
  odouble dtemp;
  gboolean wanted, doPol;
  gchar tempName[MAXKEYCHARTABLEAN+4], *ANType = "AIPS AN";
  gchar *CopyList[] = {"ARRAYX", "ARRAYY", "ARRAYZ", "GSTIA0",
		       "DEGPDY", "FREQ",   "RDATE",  "POLARX", "POLARY",
		       "UT1UTC", "DATUTC", "TIMSYS", "ARRNAM",
		       "NUMORB", "NOPCAL", "FREQID", "IATUTC", 
		       "POLTYPE", "P_REFANT", 
		       "P_DIFF01", "P_DIFF02", "P_DIFF03", "P_DIFF04", 
		       "P_DIFF05", "P_DIFF06", "P_DIFF07", "P_DIFF08", 
		       NULL};
  gchar *routine = "ObitTableANSelect";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));

  /* Poln calibration selected? */
  doPol = FALSE;
  ObitInfoListGetTest(inUV->info, "doPol", &type, (gint32*)dim, &doPol);

  /* Fully instantiate UV files */
  ObitUVFullInstantiate (inUV, TRUE, err);
  if (err->error )Obit_traceback_val (err, routine, inUV->name, retCode);
  ObitUVFullInstantiate (outUV, FALSE, err);
  if (err->error )Obit_traceback_val (err, routine, outUV->name, retCode);

  /* How many AN tables  */
  highANver = ObitTableListGetHigh (inUV->tableList, ANType);

  /* Are there any? */
  if (highANver <= 0) return OBIT_IO_OK;

  /* Loop over AN tables */
  for (iANver=1; iANver<=highANver; iANver++) {

    /* Get input table */
    numOrb  = 0;
    numPCal = 0;
    inTab = 
      newObitTableANValue (inUV->name, (ObitData*)inUV, &iANver, OBIT_IO_ReadOnly, 
			   numOrb, numPCal, err);
    if (err->error) Obit_traceback_val (err, routine, inTab->name, retCode);
    /* Find it */
    if (inTab==NULL) continue;  /* No keep looping */

    /* Open input table */
    retCode = ObitTableANOpen (inTab, OBIT_IO_ReadOnly, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inTab->name, retCode);

    /* Delete any old output table - ignore any problems */
    terr = newObitErr();
    retCode = ObitDataZapTable ((ObitData*)outUV, ANType, iANver, err);
    terr = ObitErrUnref(terr);
 
    /* Create output table */
    numOrb  = inTab->numOrb;
    numIF   = inUV->mySel->numberIF;
    if (inTab->numPCal>0) numPCal = 2 * numIF;
    else numPCal = 0;
    outTab = 
      newObitTableANValue (outUV->name, (ObitData*)outUV, &iANver, OBIT_IO_WriteOnly, 
			    numOrb, numPCal, err);
    if (err->error) Obit_traceback_val (err, routine, outUV->name, retCode);
    /* Create it? */
    Obit_retval_if_fail((outTab!=NULL), err, retCode,
			"%s: Could not create AN table %d for %s", 
			routine, iANver, outTab->name);

   /* Open output table */
  retCode = ObitTableANOpen (outTab, OBIT_IO_WriteOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inTab->name, retCode);
  
   /* Update header info */
    outTab->ArrayX  = inTab->ArrayX;
    outTab->ArrayY  = inTab->ArrayY;
    outTab->ArrayZ  = inTab->ArrayZ;
    outTab->GSTiat0 = inTab->GSTiat0;
    outTab->DegDay  = inTab->DegDay;
    outTab->Freq    = inTab->Freq;
    outTab->PolarX  = inTab->PolarX;
    outTab->PolarY  = inTab->PolarY;
    outTab->ut1Utc  = inTab->ut1Utc;
    outTab->dataUtc = inTab->dataUtc;
    outTab->FreqID  = inTab->FreqID;
    outTab->iatUtc  = inTab->iatUtc;
    outTab->P_Refant  = inTab->P_Refant;
    outTab->P_Diff01  = inTab->P_Diff01;
    outTab->P_Diff02  = inTab->P_Diff02;
    outTab->P_Diff03  = inTab->P_Diff03;
    outTab->P_Diff04  = inTab->P_Diff04;
    outTab->P_Diff05  = inTab->P_Diff05;
    outTab->P_Diff06  = inTab->P_Diff06;
    outTab->P_Diff07  = inTab->P_Diff07;
    outTab->P_Diff08  = inTab->P_Diff08;
    for (i=0; i<MAXKEYCHARTABLEAN; i++)  
      outTab->RefDate[i] = inTab->RefDate[i];
    for (i=0; i<MAXKEYCHARTABLEAN; i++)  
      outTab->TimeSys[i] = inTab->TimeSys[i];
    for (i=0; i<MAXKEYCHARTABLEAN; i++)  
      outTab->ArrName[i] = inTab->ArrName[i];
    for (i=0; i<MAXKEYCHARTABLEAN; i++)  
      outTab->polType[i] = inTab->polType[i];

   /* Copy InfoList stuff */
    ObitInfoListCopyList (inTab->myDesc->info, outTab->myDesc->info, CopyList);

    /* Poln cal info in info member */
    g_snprintf (tempName, 7, "       ");
    dim[0] = 7; dim[1] = dim[2] = dim[3] = 1;
    type = OBIT_string;
    ObitInfoListGetTest(inTab->myDesc->info, "POLTYPE", &type, dim, &tempName);
    tempName[dim[0]] = 0;
    if (doPol) {  /* If cal, blank */
      g_snprintf (tempName, 7, "       ");
    }
    ObitInfoListAlwaysPut(outTab->myDesc->info, "POLTYPE", type, dim, tempName);
    
    /* Polarization reference antenna */
    dim[0] = dim[1] = 1; type = OBIT_long; polRefAnt = 0;
    ObitInfoListGetTest(inTab->myDesc->info, "P_REFANT", &type, dim, &polRefAnt);
    ObitInfoListAlwaysPut(outTab->myDesc->info, "P_REFANT", type, dim, &polRefAnt);

    /* R-L Phase differences */
    if (!doPol) {  
      oif = 1;
      for (iif=inUV->mySel->startIF; 
	   iif<=inUV->mySel->startIF+inUV->mySel->numberIF-1;
	   iif++) {
	g_snprintf (tempName, 9, "P_DIFF%2.2d",iif);
	dim[0] = dim[1] = 1; type = OBIT_double; dtemp = 0.0;
	ObitInfoListGetTest(inTab->myDesc->info, tempName, &type, dim, &dtemp);
	g_snprintf (tempName, 9, "P_DIFF%2.2d",oif);
	ObitInfoListAlwaysPut(outTab->myDesc->info, "tempName", type, dim, &polRefAnt);
	oif++;
      }
    }

    /* Set rows */
    inRow  = newObitTableANRow (inTab);
    outRow = newObitTableANRow (outTab);
    ObitTableANSetRow (outTab, outRow, err);
    if (err->error) Obit_traceback_val (err, routine, outTab->name, retCode);

    /* Loop over table copying selected data */
    outANRow = -1;
 
    for (inANRow=1; inANRow<=inTab->myDesc->nrow; inANRow++) {
      retCode = ObitTableANReadRow (inTab, inANRow, inRow, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, inUV->name, retCode);
      if (inRow->status==-1) continue;
  
      /* Want this one? */
      wanted = ObitUVSelWantAnt(inUV->mySel, inRow->noSta);
      if (!wanted) continue;

      /* Copy selected data */
      outRow->noSta   = inRow->noSta;
      outRow->staXof  = inRow->staXof;
      outRow->PolAngA = inRow->PolAngA;
      outRow->PolAngB = inRow->PolAngB;
      for (i=0; i<8; i++) outRow->AntName[i] = inRow->AntName[i];
      for (i=0; i<3; i++) outRow->StaXYZ[i]  = inRow->StaXYZ[i];
      for (i=0; i<numOrb; i++) outRow->OrbParm[i] = inRow->OrbParm[i];
      outRow->polTypeA[0] = inRow->polTypeA[0];
      outRow->polTypeB[0] = inRow->polTypeB[0];
      /* IF dependent poln cal */
      oif = 0; 
      for (iif=inUV->mySel->startIF-1;  
	   iif<inUV->mySel->startIF+inUV->mySel->numberIF-1; 
	   iif++) { 
	if (doPol && (numPCal>0)) { /* zero cal */
	  outRow->PolCalA[2*oif]   = 0.0; 
	  outRow->PolCalA[2*oif+1] = 0.0; 
	  outRow->PolCalB[2*oif]   = 0.0; 
	  outRow->PolCalB[2*oif+1] = 0.0; 
	} else if (numPCal>0) {  /* Copy */
	  outRow->PolCalA[2*oif]   = inRow->PolCalA[2*iif]; 
	  outRow->PolCalA[2*oif+1] = inRow->PolCalA[2*iif+1];
	  outRow->PolCalB[2*oif]   = inRow->PolCalB[2*iif]; 
	  outRow->PolCalB[2*oif+1] = inRow->PolCalB[2*iif+1];
	} 
      	oif++; 
      } 

      retCode = ObitTableANWriteRow (outTab, outANRow, outRow, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, inUV->name, retCode);
    } /* end loop over rows */
    
    /* Close tables */
    retCode = ObitTableANClose (inTab, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inTab->name, retCode);
    retCode = ObitTableANClose (outTab, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, outTab->name, retCode);
 
    /* release table objects */
    inTab  = ObitTableANUnref(inTab);
    outTab = ObitTableANUnref(outTab);

    /* release row objects */
    inRow  = ObitTableANRowUnref(inRow);
    outRow = ObitTableANRowUnref(outRow);
  } /* end loop over tables */

  return retCode;
} /* end ObitTableANSelect */

