/* $Id$  */
/* Notes:
 - Not sure what the units in the Zernike polynomials are,
   for the time being use the gradients.
 - current version does reasonable things with h_ion 400 rather than
   the older 10,000 needed
 - currently antenna offset correction turned off 
 */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2013                                          */
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

#include "ObitIoN2SolNTable.h"
#include "ObitTableANUtil.h"
#include "ObitZernike.h"
#include "ObitSkyGeom.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIoN2SolNTable.c
 * ObitIoN2SolNTable function definitions.
 */

/*---------------Private function prototypes----------------*/
/** Private: Evaluate ionospheric model in direction */
static void 
ObitIoN2SolNTableModel(ObitTableNI *NITable, ObitTableNIRow *NIRow, 
		       ObitAntennaList *Ant, olong numIF, olong numPol, odouble *freqIF,
		       odouble ra, odouble dec, ofloat off[2], ofloat psign,
		       gboolean do3Dmul, gboolean doAntOff, ofloat URot3D[3][3], 
		       ObitTableSNRow *SNRow, ObitErr *err);

/** Private: Simpleminded uvw calculator */
static void
ObitIoN2SolNTableUVW (const ofloat b[3], odouble dec, ofloat ha, ofloat uvw[3]);
/*----------------------Public functions---------------------------*/

/**
 * Convert an IoN table into a SolutioN table in a given direction.
 * Output SN table will be associated with inUV.
 * Currently antenna offset form array center correction turned off 
 * \param inUV     Input uv data. 
 *                 Control info on info:
 * \li  "doAntOff" OBIT_bool (1,1,1) Correct for antenna offset from array center?
 *                 def FALSE
 * \li  "do3D"     OBIT_bool (1,1,1) Use 3D Imaging?
 *                 def TRUE
 * \param NITable  Input IonTable
 * \param outSN    if nonNULL a pointer to a previously existing Table
 * \param shift    RA and Dec offsets(deg) in which NITable to be evaluated.
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created ObitTableSN or outSN if nonNULL
 */
ObitTableSN* ObitIoN2SolNTableConvert (ObitUV *inUV, ObitTableNI *NITable,
				       ObitTableSN *outSN, ofloat shift[2], ObitErr *err)
{
  ObitIOCode retCode;
  ObitTableSNRow *SNRow   = NULL;
  ObitTableNIRow *NIRow   = NULL;
  ObitTableAN    *ANTable = NULL;
  ObitAntennaList **antennaLists=NULL;
  olong numPol, numIF, numSubA, iif, i, ierr=0;
  olong ver, numNIRow, iSNRow=0, iNIRow, iANver, highANver;
  gchar *tname;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gboolean bad, do3Dmul, allBad=TRUE, doAntOff, do3D;
  ofloat  URot3D[3][3], PRot3D[3][3], off[2], psign, fblank =  ObitMagicF();
  odouble ra, dec, raPnt, decPnt;
  gchar *routine = "ObitIoN2SolNTableConvert";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outSN;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitTableNIIsA(NITable));

  /* Antenna offset correction? */
  doAntOff = FALSE;
  ObitInfoListGetTest(inUV->info, "doAntOff", &type, dim, &doAntOff);

  /* Imaging type? */
  do3D = TRUE;
  ObitInfoListGetTest(inUV->info, "do3D", &type, dim, &do3D);

  /* create output table if needed */
  if (inUV->myDesc->jlocs>=0) numPol = MIN (2, inUV->myDesc->inaxes[inUV->myDesc->jlocs]);
  else numPol = 1;
  if (inUV->myDesc->jlocif>=0) numIF  = inUV->myDesc->inaxes[inUV->myDesc->jlocif];
  else numIF = 1;
  if (outSN==NULL) {
    tname = g_strconcat ("Calibration for: ",inUV->name, NULL);
    ver = 0;
    outSN = newObitTableSNValue (tname, (ObitData*)inUV, &ver, OBIT_IO_WriteOnly,  
				 numPol, numIF, err);
    g_free (tname);
    if (err->error) Obit_traceback_val (err, routine, inUV->name, outSN);
  }

  /* Clear any existing rows */
  ObitTableClearRows ((ObitTable*)outSN, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outSN);

  /* Get source position from header (high accuracy not needed) */
  ra  = DG2RAD * inUV->myDesc->crval[inUV->myDesc->jlocr];
  dec = DG2RAD * inUV->myDesc->crval[inUV->myDesc->jlocd];

  /* pointing position */
  raPnt  = inUV->myDesc->crval[inUV->myDesc->jlocr];
  decPnt = inUV->myDesc->crval[inUV->myDesc->jlocd];

  /* Projection onto Zernike plane */
  ObitSkyGeomRADec2Zern (raPnt, decPnt, shift[0], shift[1], &off[0], &off[1], &ierr);
  if (ierr!=0) {
    Obit_log_error(err, OBIT_Error, "%s: Error %d projecting onto Zernike Unit circle", 
		   routine, ierr);
    Obit_log_error(err, OBIT_Error, "     pos %lf %lf shift %f %f", 
		   ra, dec, shift[0], shift[1]);
    return outSN;
  }

  /* 3D/rotation matrices */
  do3Dmul = ObitUVDescShift3DPos (inUV->myDesc, shift, 0.0,  do3D, URot3D, PRot3D);

  /* Need to flip signs for GMRT & LOFAR */
  if (!strncmp(inUV->myDesc->teles, "LOFAR", 5) || 
      !strncmp(inUV->myDesc->teles, "GMRT", 4)) 
        psign = -1.0;
  else  psign =  1.0;

  /* Open IoN table for read */
  retCode = ObitTableNIOpen (NITable, OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, outSN->name, outSN);

  /* create row structure */
  NIRow = newObitTableNIRow(NITable);
  
  /* how many rows? */
  numNIRow = NITable->myDesc->nrow;

 /* Open SN table for write */
  retCode = ObitTableSNOpen (outSN, OBIT_IO_WriteOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, outSN->name, outSN);

  /* create row structure */
  SNRow = newObitTableSNRow(outSN);

   /* Attach row to output buffer */
  ObitTableSNSetRow (outSN, SNRow, err);
  if (err->error) Obit_traceback_val (err, routine, outSN->name, outSN);

  outSN->numAnt = 1;
 /* Initialize SN Row structure */
  SNRow->Time   = 0.0;
  SNRow->TimeI  = 0.0;
  SNRow->SourID = 0;
  SNRow->antNo  = 0;
  SNRow->SubA   = 0;
  SNRow->FreqID = 0;
  SNRow->IFR    = 0.0;
  SNRow->NodeNo = 1;
  SNRow->MBDelay1 = 0.0;
  SNRow->MBDelay2 = 0.0;
  SNRow->status = 1;
  for (iif = 0; iif < numIF; iif++) {
    SNRow->Real1[iif]   = 0.0;
    SNRow->Imag1[iif]   = 0.0;
    SNRow->Delay1[iif]  = 0.0;
    SNRow->Rate1[iif]   = 0.0;
    SNRow->Weight1[iif] = 0.0;
    SNRow->RefAnt1[iif] = 1;
    if (numPol>1) { /* Second polarization */
      SNRow->Real2[iif]   = 0.0;
      SNRow->Imag2[iif]   = 0.0;
      SNRow->Delay2[iif]  = 0.0;
      SNRow->Rate2[iif]   = 0.0;
      SNRow->Weight2[iif] = 0.0;
      SNRow->RefAnt2[iif] = 1;
    }
  }
  
  /* Get antenna information */
  /* How many AN tables (subarrays) */
  highANver = ObitTableListGetHigh (inUV->tableList, "AIPS AN");
  numSubA = highANver;

  /* Create AntennaLists */
  antennaLists = g_malloc0(numSubA*sizeof(ObitAntennaList*));
  for (i=0; i<numSubA; i++) antennaLists[i] = NULL;

  /* Read Info from AN tables  */
  for (i=0; i<numSubA; i++) {
    iANver = i+1;
    /* Get table */
    ANTable = newObitTableANValue (inUV->name, (ObitData*)inUV, &iANver, 
				   OBIT_IO_ReadOnly, 0, 0, 0, err);
    if ((err->error) || (ANTable==NULL)) goto cleanup;
    
    antennaLists[i] = ObitTableANGetList (ANTable, err);
    if (err->error) goto cleanup;

    /* Maximum antenna number */
    outSN->numAnt = MAX (outSN->numAnt, antennaLists[i]->number);

    /* release table object */
    ANTable = ObitTableANUnref(ANTable);
  } /* End loop over subarrays */

  /* Frequency info */
  if (inUV->myDesc->freqIF==NULL) {
    ObitUVGetFreq (inUV, err);
    if (err->error) goto cleanup;
  }

  /* Loop over NI table converting */
  for (iNIRow = 1; iNIRow<=numNIRow; iNIRow++) {
    /* read NI table */
    retCode = ObitTableNIReadRow (NITable, iNIRow, NIRow, err);
    if (err->error) goto cleanup;

    /* If entry is flagged don't bother */
    if (NIRow->status < 0) continue;

    /* Set antenna invariant values */
    SNRow->Time   = NIRow->Time;
    SNRow->TimeI  = NIRow->TimeI;
    SNRow->SourID = NIRow->SourId;
    SNRow->antNo  = SNRow->antNo;
    SNRow->SubA   = NIRow->SubA;

    /* If bad solution flag */
    bad = NIRow->weight <= 0.0;
    if (bad) { /* bad solution */
      for (iif = 0; iif < numIF; iif++) {
	SNRow->Real1[iif]   = fblank;
	SNRow->Imag1[iif]   = fblank;
	SNRow->Weight1[iif] = 0.0;
	if (numPol>1) { /* Second polarization */
	  SNRow->Real2[iif]   = fblank;
	  SNRow->Imag2[iif]   = fblank;
	  SNRow->Weight2[iif] = 0.0;
	}
      }

    } else { /* soln OK */
      for (iif = 0; iif < numIF; iif++) {
	SNRow->Weight1[iif] = 1.0;
	if (numPol>1) { /* Second polarization */
	  SNRow->Weight2[iif] = 1.0;
	}
      }
    }

    /* Are values per antenna or global? */
    if (NIRow->antNo>0) { /* one per antenna */
      SNRow->antNo  = NIRow->antNo;

      /* if this soln bad write */
      if (bad) {
	/* Write solution at beginning and end of interval */
	SNRow->Time   = NIRow->Time - 0.49 * NIRow->TimeI;
	retCode = ObitTableSNWriteRow (outSN, iSNRow, SNRow, err);
	SNRow->Time   = NIRow->Time + 0.49 * NIRow->TimeI;
	retCode = ObitTableSNWriteRow (outSN, iSNRow, SNRow, err);
      } else { /* Solution OK, calculate */

	/* convert NI model to SN */
	ObitIoN2SolNTableModel (NITable, NIRow, antennaLists[NIRow->SubA-1], 
				numIF, numPol, inUV->myDesc->freqIF, ra, dec, 
				off, psign, do3Dmul, doAntOff, URot3D, SNRow, err);

	/* write it */
 	retCode = ObitTableSNWriteRow (outSN, iSNRow, SNRow, err);
      }
      if (err->error) goto cleanup;

   } else { /* one entry for all antennas */
      /* if this soln bad write */
      if (bad) {
	/* Write solution at beginning and end of interval for each antenna */
	iSNRow = -1;
	for (i=0; i<antennaLists[NIRow->SubA-1]->number; i++) {
	  SNRow->antNo  = antennaLists[NIRow->SubA-1]->ANlist[i]->AntID;
	  SNRow->Time   = NIRow->Time - 0.49 * NIRow->TimeI;
	  retCode = ObitTableSNWriteRow (outSN, iSNRow, SNRow, err);
	  SNRow->Time   = NIRow->Time + 0.49 * NIRow->TimeI;
	  retCode = ObitTableSNWriteRow (outSN, iSNRow, SNRow, err);
	  if (err->error) goto cleanup;
	}
      } else { /* soln good - calculate antenna values */
	allBad=FALSE; /* At least this one is good */
 	for (i=0; i<antennaLists[NIRow->SubA-1]->number; i++) {
	  SNRow->antNo  = antennaLists[NIRow->SubA-1]->ANlist[i]->AntID;
	  SNRow->Time   = NIRow->Time;

 	  /* convert NI model to SN */
	  ObitIoN2SolNTableModel (NITable, NIRow, antennaLists[NIRow->SubA-1], 
				  numIF, numPol, inUV->myDesc->freqIF, ra, dec, 
				  off, psign, do3Dmul, doAntOff, URot3D,  SNRow, err);

	  /* write it */
	  iSNRow = -1;
	  retCode = ObitTableSNWriteRow (outSN, iSNRow, SNRow, err);
	  if (err->error) goto cleanup;
	}
      }
    } /* end of one entry for all antennas */

 } /* end loop over NI Table */

  /* Close tables */
  retCode = ObitTableSNClose (outSN, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;
  retCode = ObitTableNIClose (NITable, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* deallocate structures */
 cleanup:
  NIRow = ObitTableNIRowUnref(NIRow);
  SNRow = ObitTableSNRowUnref(SNRow);
  if(antennaLists) {
    for (i=0; i<numSubA; i++) 
      antennaLists[i] = ObitAntennaListUnref(antennaLists[i]);
    g_free(antennaLists);
  }
  if (err->error) Obit_traceback_val (err, routine, outSN->name, outSN);

  /* Make sure some valid solutions found */
  Obit_retval_if_fail((!allBad), err, outSN,
		      "%s: NO valid solutions found in %s", routine, NITable->name);

  return outSN;
} /* end ObitIoN2SolNTableConvert  */

/**
 * Evaluate ionspheric model in offset direction off
 * Calculate phase corrections directly from phase screen and 
 * delay correction from apparent position shift.
 * \param NITable  Input IonTable
 * \param NIRow    IonTable row
 * \param Ant      Antenna list
 * \param ra       RA of field center (rad)
 * \param dec      Dec of field center (rad)
 * \param off      RA and Dec offsets(deg/10) in which NITable to be evaluated.
 * \param psign    sign (+1.0 or -1.0) of sign and delay
 * \param do3Dmul  if 3D rotation need to be applied.
 * \param doAntOff if Antenna offset from array center correction wanted.
 * \param URot3Dl  3D rotation matrix for u,v,w 
 * \param SNRow    Output SN table row, input time and antenna used.
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created ObitTableSN or outSN if nonNULL
 */
static void 
ObitIoN2SolNTableModel(ObitTableNI *NITable, ObitTableNIRow *NIRow, 
		       ObitAntennaList *Ant, olong numIF, olong numPol, odouble *freqIF,
		       odouble ra, odouble dec, ofloat off[2], ofloat psign,
		       gboolean do3Dmul, gboolean doAntOff, ofloat URot3D[3][3], 
		       ObitTableSNRow *SNRow, ObitErr *err)
{
  /* Speed of light */
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif
  /* CI = 1/speed of light */
#ifndef CI
#define CI 1.0 / VELIGHT
#endif
/* RADSEC = earth rot rate in  rad/sec. */
#ifndef RADSEC
#define RADSEC 3.1415926535897932384 / 43200.0
#endif

  olong   j, iif;
  ofloat phase, el, az, parAng, chad, shad, ora, odec, cra, cdec;
  ofloat zaE, zaN, zaC, cE, cN, b[3], uvw[3], xyz[3];
  ofloat dx, dy;
  odouble cir, x=0.0, y=0.0, z=0.0, delayc, ratec, 
    ddec, dra, rate, pdly, dpdly, sind, cosd;
  odouble ArrLong, ArrLat, AntLst, HrAng, cosdec, sindec, darg, darg2, daz;

  cir = CI * RADSEC;

  /* get local hour angle */
  
  ArrLong = Ant->ANlist[SNRow->antNo-1]->AntLong;
  ArrLat  = Ant->ANlist[SNRow->antNo-1]->AntLat;
  AntLst = Ant->GSTIAT0 + ArrLong + NIRow->Time*Ant->RotRate;
  HrAng = AntLst - ra;

  /* assume for vla */
  chad = cos (HrAng);
  shad = sin (HrAng);

  /* Source elevation */
  cosdec = cos (dec);
  sindec = sin (dec);
  darg = sin (ArrLat) * sindec + cos (ArrLat) * cosdec * chad;
  el = (1.570796327 - acos (MIN (darg, 1.000)));

  /* Source parallactic angle */
  parAng = atan2 (cos(ArrLat) * shad, 
		  (sin(ArrLat)*cosdec - cos(ArrLat)*sindec*chad));

  /* Source azimuth */
  darg  = sindec*cos(ArrLat) - cosdec*sin(ArrLat)*chad;
  darg2 = cosdec * shad;
  daz = atan2 (darg, darg2);
  daz = fmod (daz, (2.0*G_PI));
  if (daz<0.0) daz += 2.0*G_PI;
  az = (ofloat)daz;

  /* Get zenith angle projected to east and north */
  zaE = fabs((0.5*G_PI-el)*sin(az));
  zaN = fabs((0.5*G_PI-el)*cos(az));

  /* antenna coordinates: */
  b[0] = Ant->ANlist[SNRow->antNo-1]->AntXYZ[0];
  b[1] = Ant->ANlist[SNRow->antNo-1]->AntXYZ[1];
  b[2] = Ant->ANlist[SNRow->antNo-1]->AntXYZ[2];

  /* project onto normal to source - uv aligned with East, North*/
  ObitIoN2SolNTableUVW (b, dec, HrAng, uvw);

  /* Correction for antenna offset from array center? */
  if (doAntOff) {
    /* DEBUG - force plausible height for ionosphere
       NITable->heightIon = 400.0e3;*/
    
    /* Corrected ZA to E */
    zaC = atan2(NITable->heightIon*sin(zaE)+uvw[0], NITable->heightIon);
    cE = zaC - zaE; /* angle correction to east */
    
    /* Corrected ZA to N */
    zaC = atan2(NITable->heightIon*sin(zaN)+uvw[1], NITable->heightIon);
    cN = zaC - zaN; /* angle correction to north */
    
    /* Convert antenna position offsets to screen units */
    cra  = cE * RAD2DG * 0.1;
    cdec = cN * RAD2DG * 0.1;
  } else { /* No antenna offset corrections */
    cra = 0.0;
    cdec = 0.0;
  } /* End antenna offset correction */

  /* delay and rate in sec and  sec/sec. (want corrections). 
     the formulae are written  for the right hand coordinate  system. */
  /* uncorrected rate */
  x = b[0]; y = b[1]; z = b[2];
  /*DEBUG delay = CI * ( (x * chad - y * shad) * cosdec + z * sindec);*/
  rate  = cir * (-x * shad - y * chad) * cosdec;

  /* position shift in this direction  */
  dra  =  0.0;
  ddec =  0.0;

  /* evaluate zernike model for apparent position shift */
  ora  = off[0] + cra;
  odec = off[1] + cdec;
  for (j=0; j<NITable->numCoef; j++) {  /* loop  L100: */
    dra  += NIRow->coef[j]*ObitZernikeGradX(j+2, ora, odec);
    ddec += NIRow->coef[j]*ObitZernikeGradY(j+2, ora, odec);
  }  /* end loop  L100: */
  
  /* corrected shift parameters with 3D shift */
  dx = -dra  * DG2RAD;
  dy = -ddec * DG2RAD;
  if (do3Dmul) {
    xyz[0] = dx*URot3D[0][0] + dy*URot3D[1][0];
    xyz[1] = dx*URot3D[0][1] + dy*URot3D[1][1];
    xyz[2] = dx*URot3D[0][2] + dy*URot3D[1][2];
  } else {
    xyz[0] = dx;
    xyz[1] = dy;
    xyz[2] = 0.0;
  }

  /* source position error as corrections in radians */
  dra  = -dra  * DG2RAD;
  ddec = -ddec * DG2RAD; 

  /* correct hour angle */
  chad = cos (HrAng - dra);
  shad = sin (HrAng - dra);

  /* correct declination */
  sind = sin (dec + ddec);
  cosd = cos (dec + ddec);

  /* corrected delay and rate - 
     for delay only need component due to apparent position shift */
  /*DEBUGdelayc = CI * ( (x * chad - y * shad) * cosd + z * sind);*/
  delayc = CI*(uvw[0]*xyz[0] + uvw[1]*xyz[1] + uvw[2]*xyz[2]);
  ratec  = cir * (-x * shad - y * chad) * cosd;

  /* correction of delay (turns/Hz) and rate */
  pdly  = delayc;
  dpdly = ratec - rate;

  /* set corrections */
  for (iif = 0; iif < numIF; iif++) {
    phase = pdly * 2.0 * G_PI * freqIF[iif];
    SNRow->Real1[iif]   = cos (phase);
    SNRow->Imag1[iif]   = psign*sin (phase);
    SNRow->Delay1[iif]  = psign*pdly;
    SNRow->Rate1[iif]   = 0.0;
    if (numPol>1) { /* Second polarization */
      SNRow->Real2[iif]   = cos (phase);
      SNRow->Imag2[iif]   = psign*sin (phase);
      SNRow->Delay2[iif]  = psign*pdly;
      SNRow->Rate2[iif]   = 0.0;
    }
  }

} /* end ObitIoN2SolNTableModel */

/**
 * UVW for given baseline and geometry
 * \param b   Baseline vector, (VLA Convention)
 * \param dec Source declination (radians)
 * \param ha  Source hour angle (radians)
 * \param uvw U, V, W in same units as B
 */
static void
ObitIoN2SolNTableUVW (const ofloat b[3], odouble dec, ofloat ha, ofloat uvw[3])
{
  odouble cosdec, sindec;
  ofloat     sinha, cosha, vw;

  cosdec = cos (dec);
  sindec = sin (dec);
  cosha = cos (ha);
  sinha = sin (ha);
  vw = b[0]*cosha - b[1]*sinha;
  uvw[0] =  b[0]*sinha + b[1]*cosha;
  uvw[1] = -vw*sindec + b[2]*cosdec;
  uvw[2] =  vw*cosdec + b[2]*sindec;
}
/* end ObitIoN2SolNTableUVW */
