/* $Id$  */
/* this version        2008-10-01 20:20:00  juan.uson      */
/* J1 extended with large angle approximation              */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2012                                          */
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

#include "ObitPBUtil.h"
#include "ObitSkyModel.h"
#include <math.h>

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPBUtil.c
 * ObitPBUtil function definitions.
 *
 * Antenna primary beam shape utility class.
 */

/** AIPSish Primary beam calculation */
ofloat pbfact (olong pbtype, olong nfreq, odouble *pbfreq, ofloat pbfsiz, 
	       olong ifreq, ofloat radius);
/*----------------------Public functions---------------------------*/

/**
 * Compute VLA beam shape from a fitted polynomial
 * From the AIPSish $APLSUB/PBCALC.FOR
 * \param Angle  Angle from the pointing position (deg)
 * \param Freq   Frequency (Hz) of observations
 * \param pbmin  Minimum antenna gain 0=>0.01
 * \return Fractional antenna power [pbmin, 1]
 */
ofloat ObitPBUtilPoly (odouble Angle, odouble Freq, ofloat pbmin)
{
  ofloat bmfact;
  olong  i;
  odouble x, bm[7], bmult;
  static ofloat  table[8][3] = {
    /* {-0.897e-3,  2.71e-7 , -0.242e-10}, Rick's value */
    {-1.051e-2,  4.276e-7, -5.380e-11}, /* Fitted to VLSS fields */
    {-0.935e-3,  3.23e-7 , -0.378e-10},
    {-1.343e-3,  6.579e-7, -1.186e-10},
    {-1.372e-3,  6.940e-7, -1.309e-10},
    {-1.306e-3,  6.253e-7, -1.100e-10},
    {-1.305e-3,  6.155e-7, -1.030e-10},
    {-1.417e-3,  7.332e-7, -1.352e-10},
    {-1.321e-3,  6.185e-7, -0.983e-10}};

  /* Hack - use old AIPSish routine 
  bmfact = pbfact (2, 1, &Freq, 24.5, 0, (float)Angle);
  return bmfact;*/
  
  bmfact = 0.0;
  /* which VLA band */
  bmult = Freq * 60.0e-9;
  if (Freq < 0.15) {
    i = 0;
  } else if (Freq < 1.1e9) {
    i = 1;
  } else if (Freq < 2.0e9) {
    i = 2;
  } else if (Freq < 6.0e9) {
    i = 3;
  } else if (Freq < 10.0e9) {
    i = 4;
  } else if (Freq < 18.0e9) {
    i = 5;
  } else if (Freq < 30.0e9) {
    i = 6;
  } else {
    i = 7;
  } 
  bm[0] = table[i][0];
  bm[1] = table[i][1];
  bm[2] = table[i][2];
  bm[3] = 0.0e0;
  bm[4] = 0.0e0;

  x = (Angle * bmult) * (Angle * bmult);

  bmfact = 1.0 + (bm[0] + (bm[1] + (bm[2] + (bm[3] + bm[4] * x) * x) * x) * x) * x;
  
  /* here, make some "reasonable" estimate on the rumbling around 
     in far sidelobes, and the depths of the nulls... */
  if (pbmin<=0.0) pbmin = 0.01;
  bmfact = MIN (1.0, MAX (bmfact, pbmin));

  return bmfact;
} /*  end ObitPBUtilPoly */

/**
 * Compute Antenna beam shape assuming uniform illumination of an antenna
 * with diameter antSize.  The power pattern is calculated 
 * from the pointing position and for observing frequency freq (Hz). 
 * The power pattern (2 * j1(x) / x) ** 2 of a uniformly illuminated 
 * circular aperture is used, since it fits the observations better 
 * than the standard PBCOR beam does.  If the relative gain is less 
 * than pbmin, it is set to pbmin. 
 * vscale is a measured constant inversely proportional to the 
 * VLA primary beamwidth, which is assumed to scale as 1./freq. 
 * vscale = 4.487e-9 corresponds to a 29.4 arcmin fwhm at 1.47 ghz. 
 * the actual scale is determined from the antenna size (antSize). 
 * xmax = value of x where the series approximation to the J1 goes from the
 * small angle approximation to the large angle approximation.
 * Note: this routine is probably only useful for the VLA but might 
 * be ok for a homogenous array of uniformly illuminated antennas where 
 * the beam scales from the VLA beam by the ratio of antenna diameters. 
 * From the AIPSish $FOURMASS/SUB/PBUTIL.FOR
 * \param Angle   Angle from the pointing position (deg)
 * \param Freq    Frequency (Hz) of observations
 * \param antSize Antenna diameter in meters. (defaults to 25.0)
 * \param pbmin   Minimum antenna gain 0=>0.05
 * \return Fractional antenna power [pbmin, 1]
 */
ofloat ObitPBUtilJinc (odouble Angle, odouble Freq, ofloat antSize, 
		       ofloat pbmin)
{
  ofloat bmfact = 1.0;
  ofloat  x, u, scale, pb, asize, xx, f, t;
  /* coefficients c from Abramowitz and Stegun, eq. 9.4.4 */
  static ofloat c1 = -0.56249985;
  static ofloat c2 =  0.21093573;
  static ofloat c3 = -0.03954289;
  static ofloat c4 =  0.00443319;
  static ofloat c5 = -0.00031761;
  static ofloat c6 =  0.00001109;
  /* coefficients d and e from Abramowitz and Stegun, eq. 9.4.6 */
  
  static ofloat d1 =  0.79788456;
  static ofloat d2 =  0.00000156;
  static ofloat d3 =  0.01659667;
  static ofloat d4 =  0.00017105;
  static ofloat d5 = -0.00249511;
  static ofloat d6 =  0.00113653;
  static ofloat d7 = -0.00020033;

  static ofloat e1 = -2.35619449;
  static ofloat e2 =  0.12499612;
  static ofloat e3 =  0.00005650;
  static ofloat e4 = -0.00637879;
  static ofloat e5 =  0.00074348;
  static ofloat e6 =  0.00079824;
  static ofloat e7 = -0.00029166;

  static ofloat vscale=  4.487e-9;
  
  /* transition to large angle J1 */
  static ofloat xmax = 3.0;

  /* default antenna size */
  asize = antSize;
  if (asize <= 0.0) asize = 25.0;

  /* beam scale size at 1.47 GHz */
  scale = vscale * 25.0 / asize;

  x = scale * Angle * Freq * DG2RAD;
  if (x  <  xmax) {
    u = x * x / 9.0;
    pb = 0.5 + u*(c1 + u*(c2 + u*(c3 + u*(c4 + u*(c5 + u*c6)))));
  } else {
    x = MIN (x, 5*xmax);
    xx = x / 3.0;
    u = 3.0 / x;
    f = d1 + xx*(d2 + xx*(d3 + xx*(d4 + xx*(d5 + xx*(d6 + xx*d7)))));
    t = x + e1 + u*(e2 + u*(e3 + u*(e4 + u*(e5 + u*(e6 + u*e7)))));
    pb = f * cos(t) / (pow(x,1.5));
  }	

  /* Allow going to 0.001 if using J1 beam */
  bmfact = MAX (0.001, 4.* pb * pb);

  /* here, make some "reasonable" estimate on the rumbling around 
     in far sidelobes, and the depths of the nulls... */
  if (pbmin<=0.0) pbmin = 0.01;
  bmfact = MIN (1.0, MAX (bmfact, pbmin));

  return bmfact;
} /*  end ObitPBUtilJinc */

/**
 * Calculates the relative gain at a reference frequency (refFreq) 
 * relative to the average of a set of frequencies (Freq) for a given
 * offset from the antenna pointing position (Angle).
 * Uses ObitPBUtilPoly (VLA assumed) for frequencies < 1 GHz
 * and ObitPBUtilJinc at higher frequencies
 * Adopted from the AIPSish $FOURMASS/SUB/PBUTIL.FOR PBFACT
 * \param Angle   Angle from the pointing position (deg)
 * \param nfreq   number of frequencies in Freq
 * \param Freq    Frequencies (Hz) of observations.
 * \param antSize Antenna diameter in meters. (defaults to 25.0)
 * \param pbmin   Minimum antenna gain Jinc 0=>0.05, poly 0=> 0.01
 * \param refFreq Reference frequency (Hz) for which rel. gain is desired
 * \return Relative gain at freq refFreq wrt average of Freq.
 */
ofloat ObitPBUtilRelPB (odouble Angle, olong nfreq, odouble *Freq, ofloat antSize, 
			ofloat pbmin, odouble refFreq)
{
  ofloat PBfact;
  olong i;
  ofloat PBref, iPBref, sum, pb;
  gboolean doJinc;
  
  /* Which beam shape function to use? */
  doJinc = (Freq[0] >= 1.0e9);
  
  /* Gain at refFreq */
  if (doJinc) PBref = ObitPBUtilJinc(Angle, refFreq, antSize, pbmin);
  else        PBref = ObitPBUtilPoly(Angle, refFreq, pbmin);
  
  /* inverse of PBref */
  if (PBref>0.0) iPBref = 1.0 / PBref;
  else iPBref = 1.0;
  
  sum = 0.0;
  /* sum gain relative to ref. freq. */
  for (i=0; i<nfreq; i++) {
    if (doJinc) pb = ObitPBUtilJinc(Angle, Freq[i], antSize, pbmin);
    else        pb = ObitPBUtilPoly(Angle, Freq[i], pbmin);
    sum += pb * iPBref;  /* Sum of rel frequency gain */
  }
  /*  Compute average relative gain */
  if (sum<=0.0)  PBfact = 1.0;
  else  PBfact = nfreq  / sum;
  
  return PBfact;
} /* end ObitPBUtilRelPB */

/**
 * Calculates the correction needed to undo the effect of a 
 * on amplitude.  Note: this is the power pattern.
 * Uses ObitPBUtilPoly (VLA assumed) for frequencies < 1 GHz
 * and ObitPBUtilJinc at higher frequencies
 * Adopted from the AIPSish $FOURMASS/SUB/PBUTIL.FOR PBFACT
 * \param Angle   Intended angle from the intended pointing position (deg)
 * \param AngleO  Actual angle from the intended pointing position (deg)
 * \param antSize Antenna diameter in meters. (defaults to 25.0)
 * \param pbmin   Minimum antenna gain Jinc 0=>0.05, poly 0=> 0.01
 * \param Freq    Frequency (Hz) for which rel. gain is desired
 * \return amplitude correction
 */
ofloat ObitPBUtilPntErr (odouble Angle, odouble AngleO, ofloat antSize, 
			 ofloat pbmin, odouble Freq)
{
  ofloat PBfact;
  ofloat PBref, PBoff;
  gboolean doJinc;
  
  /* Which beam shape function to use? */
  doJinc = (Freq >= 1.0e9);
  
  /* Gain at correct position */
  if (doJinc) PBref = ObitPBUtilJinc(Angle, Freq, antSize, pbmin);
  else        PBref = ObitPBUtilPoly(Angle, Freq, pbmin);
  
   /* Gain at offset position */
  if (doJinc) PBoff = ObitPBUtilJinc(AngleO, Freq, antSize, pbmin);
  else        PBoff = ObitPBUtilPoly(AngleO, Freq, pbmin);

  PBfact = PBref/PBoff;
  
  return PBfact;
} /* end ObitPBUtilPntErr */

/**
 * Derive an ObitTableCC from the input one in which the fluxes
 * are corrected by the relative antenna gains between refFreq and 
 * the average of Freq.
 * Also processes CC tables with tabulated spectra, i.e.
 * Param[3] between 20 and 29, each channel is processed.
 * From the AIPSish $FOURMASS/SUB/PBUTIL.FOR PBFCCT 
 * \param image    input image with input CC table
 * \param inCCver  input CC table
 * \param outCCver Desired output CC table on image, if 0 then new
 *                 value used returned.
 * \param nfreq    number of frequencies in Freq
 * \param Freq     Frequencies (Hz) of observations.
 * \param antSize  Antenna diameter in meters. (defaults to 25.0)
 * \param pbmin   Minimum antenna gain Jinc 0=>0.05, poly 0=> 0.01
 * \param refFreq  Reference frequency (Hz) for which CC table is needed
 * \param startCC  [in] the desired first CC number (1-rel)
 *                 [out] the actual first CC number in returned table
 * \param endCC    [in] the desired highest CC number, 0=> to end of table
 *                 [out] the actual highest CC bumber in returned table
 * \param err      Obit error/message object
 * \return pointer to ObitCCtable, Unref when no longer needed.
 */
ObitTableCC *ObitPBUtilCCCor(ObitImage *image, olong inCCver, olong *outCCver, 
			     olong nfreq, odouble *Freq, ofloat antSize, ofloat pbmin,
			     odouble refFreq,  olong *startCC, olong *endCC,
			     ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTable *tempTable=NULL;
  ObitTableCC *inCCTable = NULL, *outCCTable = NULL;
  ObitTableCCRow *CCRow = NULL;
  gchar *tabType = "AIPS CC";
  odouble *specFreq=NULL;
  ofloat Angle;
  olong i, j, ver, irow, orow, nSpec, offset=0;
  gchar keyword[20];
  gboolean isSpec=FALSE;
  ObitSkyModelCompType modType;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "ObitPBUtilCCCor";

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outCCTable;
  g_assert (ObitImageIsA(image));

  /* Get CC table */
  ver = inCCver;
  tempTable = newObitImageTable (image,OBIT_IO_ReadOnly, tabType, &ver, err);
  if ((tempTable==NULL) || (err->error)) 
     Obit_traceback_val (err, routine, image->name, outCCTable);
  inCCTable = ObitTableCCConvert(tempTable);
  tempTable = ObitTableUnref(tempTable);
  if (err->error) Obit_traceback_val (err, routine, image->name, outCCTable);
  
  /* Open input */
  retCode = ObitTableCCOpen (inCCTable, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, image->name, outCCTable);
  /* Create table row */
  CCRow = newObitTableCCRow (inCCTable);
  isSpec = inCCTable->noParms>4;  /* Might this have a tabulated spectrum? */
  
  /* Create output CC table */
  ver = *outCCver;
  outCCTable = newObitTableCCValue ("PB Corrected", (ObitData*)image,
				    &ver, OBIT_IO_WriteOnly, inCCTable->noParms, 
				    err);
  if (err->error) Obit_traceback_val (err, routine, image->name, outCCTable);
  *outCCver = ver;  /* save if defaulted (0) */
      
  /* Open output */
  retCode = ObitTableCCOpen (outCCTable, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, image->name, outCCTable);
  
  /* loop over table */
  orow = 0;
  for (j=*startCC; j<=*endCC; j++) {
    irow = j;
    retCode = ObitTableCCReadRow (inCCTable, irow, CCRow, err);
    if (retCode == OBIT_IO_EOF) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: EOF CC table row %d in %s", 
		     routine, irow, image->name);
    }
    if  (err->error) Obit_traceback_val (err, routine, image->name, outCCTable);

    /* Determine angle wrt pointing position */
    Angle = ObitImageDescAngle(image->myDesc, CCRow->DeltaX, CCRow->DeltaY);

    /* Correct flux density */
    CCRow->Flux *= ObitPBUtilRelPB ((odouble)Angle, nfreq, Freq, antSize, refFreq, 
				    pbmin);

    /* Is this a tabulated spectrum CC? */
    if (isSpec) {
      if ((inCCTable->noParms>4) && (CCRow->parms[3]>=19.99) && 
	  (CCRow->parms[3]<=29.99)) {
	/* Need to initialize? */
	if ((nSpec<=0) || (specFreq==NULL)) {
	  nSpec = 0;
	  ObitInfoListGetTest(image->myDesc->info, "NSPEC", &type, dim, &nSpec);
	  /* get number of and channel frequencies for CC spectra from 
	     CC table on first image in mosaic */
	  if (nSpec>0) {
	    specFreq = g_malloc0(nSpec*sizeof(odouble));
	    for (i=0; i<nSpec; i++) {
	      specFreq[i] = 1.0;
	      sprintf (keyword, "FREQ%4.4d",i+1);
	      ObitInfoListGetTest(image->myDesc->info, keyword, &type, dim, &specFreq[i]);
	    } /* end loop reading frequencies */
	    /* Get model type and offset in record of start of spectrum */
	    modType = CCRow->parms[3] + 0.5;
	    /* Offset in record */
	    offset = 4;
	    /* end initialize */
	  } else {isSpec = FALSE; nSpec = 0;}
	} /* end initialize */

	/* Correct tabulated spectrum */
	for (i=0; i<nSpec; i++) {
	  CCRow->parms[offset+i] *=  ObitPBUtilRelPB ((odouble)Angle, 1, &specFreq[i], 
						     antSize, refFreq, pbmin);
	}
      } else isSpec = FALSE;  /* No tabulated spectrum */
    } /* end tabulated spectrum  */

    /* Write output */
    orow++;
    retCode = ObitTableCCWriteRow (outCCTable, orow, CCRow, err);
    if  (err->error) Obit_traceback_val (err, routine, image->name, outCCTable);
  } /* end loop over table */

  /* Close */
  retCode = ObitTableCCClose (outCCTable, err);
  retCode = ObitTableCCClose (inCCTable, err);
  if  (err->error) Obit_traceback_val (err, routine, image->name, outCCTable);
  inCCTable = ObitTableUnref(inCCTable);
  CCRow = ObitTableRowUnref(CCRow);
  if (specFreq) g_free(specFreq);

  /* Set actual values in output table */
  *startCC = 1;
  *endCC   = orow;

 return outCCTable;
} /* end ObitPBUtilCCCor */

/**
 * Derive an image (FArray) from the input one in which the pixels
 * are corrected by the relative antenna gains between refFreq and 
 * the average of Freq.
 * From the AIPSish $FOURMASS/SUB/PBUTIL.FOR PBFSCI
 * \param inImage  input image with
 * \param inPlane  Desired plane in inImage, 1-rel pixel numbers on planes 3-7; 
 *                 ignored if memOnly
 * \param nfreq    number of frequencies in Freq
 * \param Freq     Frequencies (Hz) of observations.
 * \param antSize  Antenna diameter in meters. (defaults to 25.0)
 * \param pbmin    Minimum antenna gain Jinc 0=>0.05, poly 0=> 0.01
 * \param refFreq  Reference frequency (Hz) for which CC table is needed
 * \param err      Obit error/message object
 * \return pointer to ObitFArray, Unref when no longer needed.  NULL on error
 */
ObitFArray* ObitPBUtilImageCor(ObitImage *inImage, olong *inPlane, 
			       olong nfreq, odouble *Freq, 
			       ofloat antSize, ofloat pbmin, odouble refFreq, 
			       ObitErr *err)
{
  ObitIOSize IOBy;
  olong blc[IM_MAXDIM], trc[IM_MAXDIM];
  gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong ix, iy, indx, pos[2];
  ofloat inPixel[2], *out;
  ofloat DeltaX, DeltaY, Angle, fblank = ObitMagicF();
  ObitImageDesc *inDesc;
  ObitFArray *outFA=NULL;
  gchar *routine = "ObitPBUtilImageCor";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outFA;
  g_assert (ObitImageIsA(inImage));

  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;

  /* Do I/O by plane and all of plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListPut (inImage->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
  dim[0] = 7;
  for (i=0; i<5; i++) blc[i+2] = trc[i+2] = inPlane[i];
  ObitInfoListPut (inImage->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (inImage->info, "TRC", OBIT_long, dim, trc, err);
  if (err->error) Obit_traceback_val (err, routine, inImage->name, outFA);

  /* Open image */
  ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_val (err, routine, inImage->name, outFA);

  /* Read input plane */
  ObitImageRead (inImage, NULL , err);
  if (err->error) Obit_traceback_val (err, routine, inImage->name, outFA);

  /* Get aray pointer - use this array as return value */
  outFA = ObitFArrayRef(inImage->image);
  pos[0] = pos[1] = 0;
  out = ObitFArrayIndex (outFA, pos);

  inDesc = inImage->myDesc; /* Input descriptor */

  /* Loop over image  */
  for (iy = 1; iy<=inDesc->inaxes[1]; iy++) { /* loop in y */
    inPixel[1] = (ofloat)iy;
    /* Get offset from reference position */
    DeltaY = (inPixel[1] - inDesc->crpix[1]) * inDesc->cdelt[1];
    for (ix = 1; ix<=inDesc->inaxes[0]; ix++) {/* loop in x */
      inPixel[0] = (ofloat)ix;

      /* array index in out for this pixel */
      indx = (iy-1) * inDesc->inaxes[0] + (ix-1);

      /* Is this pixel valid? */
      if (out[indx] != fblank) {

	/* Get offset from reference position */
	DeltaX = (inPixel[0] - inDesc->crpix[0]) * inDesc->cdelt[0];

	/* Determine angle wrt pointing position */
	Angle = ObitImageDescAngle(inDesc, DeltaX, DeltaY);

	/* Make correction */
	out[indx] *= ObitPBUtilRelPB ((odouble)Angle, nfreq, Freq, antSize, refFreq, pbmin); 
      } /* end if pixel valid */

    } /* end loop over x */
  } /* end loop over y */
  

  /* Close */
  ObitImageClose (inImage, err);
  if (err->error) Obit_traceback_val (err, routine, inImage->name, outFA);

  /* Detatch FArray from image */
  inImage->image = ObitFArrayUnref(inImage->image);

  return outFA;
} /* end ObitPBUtilImageCor */

/**
 * Calculates the relative gain (normalized to unity at the pointing
 * position) of the primary beam at angular offset RADIUS (deg)
 * from the pointing position and for observing frequency FREQ (Hz).
 * The power pattern (2 * J1(X) / X) ** 2 of a uniformly illuminated 
 * circular aperture is used, since it fits the observations better
 * than the standard PBCOR beam does.  If the relative gain is less
 * than pbmin = 0.05, it is set to pbmin.
 *    vscale is a measured constant inversely proportional to the
 * VLA primary beamwidth, which is assumed to scale as 1./freq.
 * vscale = 4.487E-9 corresponds to a 29.4 arcmin fwhm at 1.47 GHz.
 * The actual scale is determined from the antenna size (pbfsiz).
 * xmax = value of x yielding pb = pbmin = 0.05, beyond which the 
 * series approximation loses accuracy.
 *    NOTE: This routine is probably only useful for the VLA but might
 * be OK for a homogenous array of uniformly illuminated antennas where
 * the beam scales from the VLA beam by the ratio of antenna diameters.
 * Translated from the AIPSISH PBUTIL.FOR:PBFACT
 * \param  pbtype    Primary Beam correction type 1= rel(default), 2=abs
 * \param  nfreq     Number of frequencies in pbfreq 
 * \param  pbfreq    Frequencies (Hz) going into the average.
 * \param  pbfsiz    Antenna diameter (m)
 * \param  ifreq     0-rel Index in pbfreq of current frequency
 * \param  radius    Distance from pointing center in deg.
 * \return  Relative or absolute gain
 */
ofloat pbfact (olong pbtype, olong nfreq, odouble *pbfreq, ofloat pbfsiz, 
	       olong ifreq, ofloat radius) {
  ofloat pbf = 1.0;
  olong   i;
  ofloat      *p=NULL, sum, x, u, scale, pb, asize;
  /* Coefficients C from Abramowitz  and Stegun, eq. 9.4.4 */
  static odouble c1=-0.56249985, c2=0.21093573, c3=-0.03954289;
  static odouble c4=0.00443319, c5=-0.00031761, c6=0.00001109;
  static ofloat pbmin=0.05, vscale=4.487e-9, xmax=3.00751;

  asize = pbfsiz;
  if (asize <= 0.0) asize = 25.0;
  /* Beam scale size at 1.47 GHz */
  scale = vscale * asize / 25.0;
  sum = 0.0;
  if (pbtype == 2) {
    /* Absolute gain */
    x = scale * radius * pbfreq[ifreq];
    if (x  <  xmax) {
      u = x * x / 9.;
      pb = 0.5 + u*(c1 + u*(c2 + u*(c3 + u*(c4 + u*(c5 + u*c6)))));
      pbf = 4.* pb * pb;
    } else {
      pbf = pbmin;
    } 
  } else {
    /* Gain relative to avg. freq 
       Compute antenna power gains */
    p = g_malloc0(nfreq*sizeof(ofloat));
    for (i=0; i<nfreq; i++) {
      x = scale * radius * pbfreq[i];
      if (x  <  xmax) {
	u = x * x / 9.;
	pb = 0.5 + u*(c1 + u*(c2 + u*(c3 + u*(c4 + u*(c5 + u*c6)))));
	pb = 4.* pb * pb;
      } else {
	pb = pbmin;
      } 
      p[i] = pb;
      sum = sum + p[i-1];
    } /* end loop over frequency */
    /* Compute relative gain */
    if (sum <= 0.0) {
      pbf = 1.0;
    } else {
      pbf = nfreq * p[ifreq] / sum;
    } 
    if (p) g_free(p);
  } /* end rel gain */
  return pbf;
} /* end of routine pbfact */ 


