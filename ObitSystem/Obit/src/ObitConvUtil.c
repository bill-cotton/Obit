/* $Id: ObitConvUtil.c,v 1.2 2007/08/31 17:24:03 bcotton Exp $ */
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

#include "ObitConvUtil.h"
#include "ObitFeatherUtil.h"
#include "ObitImageUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitConvUtil.c
 * ObitConvUtil utility function definitions.
 *
 * Routines for image convolution
 */


/*----------------------Private functions---------------------------*/
/** Private: */

/*----------------------Public functions---------------------------*/

/**
 * (de)Convolve an Image with an FArray and write outImage  
 *  This routine convolves all selected planes in inImage with convFn if 
 *  doDivide is FALSE, else it does a linear deconvolution
 *  Operations are performed using FFTs
 * \param inImage   Input ObitImage 
 * \param convFn    Convolving Function    
 * \param doDovide  If true divide FT of convFn into FT of inImage,
 *                  else multiply
 * \param rescale   multiplication factor to scale output to correct units
 * \param outImage  Output ObitImage must be a clone of inImage
 *                  Actual convolution size must be set externally 
 * \param err       ObitErr for reporting errors.
 */
void ObitConvUtilConv (ObitImage *inImage, ObitFArray *convFn, 
		       gboolean doDivide, ofloat rescale,
		       ObitImage *outImage, ObitErr *err)
{
  ObitIOCode   iretCode, oretCode;
  olong      ndim=2, naxis[2], blc[2], trc[2], cen[2];
  ofloat Beam[3];
  ObitFFT    *FFTfor=NULL, *FFTrev=NULL;
  ObitFArray *padConvFn=NULL, *padImage=NULL, *tmpArray=NULL;
  ObitCArray *wtArray=NULL, *FTArray=NULL;
  gchar *routine = "ObitConvUtilConv";

  if (err->error) return;  /* existing error? */

  /* Create FFTs */
  FFTfor = ObitFeatherUtilCreateFFT(inImage, OBIT_FFT_Forward);
  FFTrev = ObitFeatherUtilCreateFFT(inImage, OBIT_FFT_Reverse);

  /* Pad convolving function */
  naxis[0] = FFTfor->dim[0];  naxis[1] = FFTfor->dim[1]; 
  padConvFn = ObitFArrayCreate("Pad Conv Fn", ndim, naxis);
  ObitFeatherUtilPadArray (FFTfor, convFn, padConvFn);

  /* FFT Convolving function to wtArray */
  wtArray = ObitFeatherUtilCreateFFTArray (FFTfor);
  ObitFArray2DCenter (padConvFn); /* Swaparoonie to FFT order */
  ObitFFTR2C (FFTfor, padConvFn, wtArray);

  /* Pad array for image */
  naxis[0] = FFTfor->dim[0];  naxis[1] = FFTfor->dim[1]; 
  padImage = ObitFArrayCreate("Pad Image", ndim, naxis);
  FTArray  = ObitFeatherUtilCreateFFTArray (FFTfor);

  /* Open input image */
  iretCode = ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  if (err->error) goto cleanup;

  /* Save output resolution */
  Beam[0] = outImage->myDesc->beamMaj;
  Beam[1] = outImage->myDesc->beamMin;
  Beam[2] = outImage->myDesc->beamPA;

  /* Copy descriptor */
  outImage->myDesc = ObitImageDescCopy(inImage->myDesc, outImage->myDesc, err);
  if (err->error) goto cleanup;

  /* Set output resolution */
  outImage->myDesc->beamMaj = Beam[0];
  outImage->myDesc->beamMin = Beam[1];
  outImage->myDesc->beamPA  = Beam[2];

  /* Force to float pixels */
  outImage->myDesc->bitpix=-32;
  outImage->myDesc->maxval = -1.0e20;
  outImage->myDesc->minval =  1.0e20;

  /* Open output image - Use external buffer for writing output */
  outImage->extBuffer = TRUE;
  oretCode = ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  if (err->error) goto cleanup;

  /* Normalize rescale for FFT */
  rescale /= (ofloat)(FFTfor->dim[0] * FFTfor->dim[1]);

  /* Loop over planes until hitting EOF */
  while (iretCode!= OBIT_IO_EOF) {
    iretCode = ObitImageRead (inImage, NULL, err);
    if (iretCode == OBIT_IO_EOF) break;
    if (err->error) goto cleanup;

    /* Pad image */
    ObitFeatherUtilPadArray (FFTfor, inImage->image, padImage);
    
    /* FFT Convolving function to FTArray */
    ObitFArray2DCenter (padImage); /* Swaparoonie to FFT order */
    ObitFFTR2C (FFTfor, padImage, FTArray);

    /* Multiply or divide by transfer function */
    if (doDivide) ObitCArrayDiv (FTArray, wtArray, FTArray);
    else ObitCArrayMul (FTArray, wtArray, FTArray);

    /* Back FFT */
    ObitFFTC2R(FFTrev, FTArray, padImage);
    ObitFArray2DCenter (padImage);/* Swaparoonie */

    /* DEBUG 
    ObitImageUtilArray2Image ("ConvolDebug1.fits",1,padImage, err);*/

    /* Get window to extract */
    cen[0] = FFTfor->dim[0]/2;  cen[1] = FFTfor->dim[1]/2; 
    blc[0] = cen[0] - inImage->image->naxis[0] / 2; 
    /*trc[0] = cen[0] - 1 + inImage->image->naxis[0] / 2;*/
    trc[0] = cen[0] + inImage->image->naxis[0] / 2;
    trc[0] -= (trc[0]-blc[0]+1) - inImage->image->naxis[0];
    blc[1] = cen[1] - inImage->image->naxis[1] / 2; 
    /*trc[1] = cen[1] - 1 + inImage->image->naxis[1] / 2;*/
    trc[1] = cen[1] + inImage->image->naxis[1] / 2;
    trc[1] -= (trc[1]-blc[1]+1) - inImage->image->naxis[1];
    
    /* Extract */
    tmpArray = ObitFArraySubArr(padImage, blc, trc, err);
    if (err->error) goto cleanup;

    /* rescale units */
    ObitFArraySMul (tmpArray, rescale);

    /* DEBUG
    ObitImageUtilArray2Image ("ConvolDebug1.fits",1,tmpArray, err); */

    /* Write plane */
    oretCode = ObitImageWrite(outImage, tmpArray->array, err);
    if (err->error) goto cleanup;
    tmpArray  = ObitFArrayUnref(tmpArray);
  } /* end loop over input planes */


  /* Close input */
  iretCode = ObitImageClose (inImage, err);
  if (err->error) goto cleanup;
  /* Free image buffer */
  inImage->image = ObitFArrayUnref(inImage->image);

  /* Close output */
  oretCode = ObitImageClose (outImage, err);
  if (err->error) goto cleanup;
  /* Unset external buffer for writing */
  outImage->extBuffer = FALSE;

 cleanup:
  wtArray   = ObitCArrayUnref(wtArray);
  FTArray   = ObitCArrayUnref(FTArray);
  padConvFn = ObitFArrayUnref(padConvFn);
  padImage  = ObitFArrayUnref(padImage);
  tmpArray  = ObitFArrayUnref(tmpArray);
  FFTfor    = ObitFFTUnref(FFTfor);
  FFTrev    = ObitFFTUnref(FFTrev);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
} /* end ObitConvUtilConv */

/**
 * Create an ObitFArray containing a unit area Gaussian in the center
 * \param inImage ObitImage giving the geometry of the output array
 * \param Beam    Gaussian major, minor axis, position angle all in deg.
 *                If there is no pixel spacing information in the inImage
 *                descriptor then the size is assumed in pixels.
 * \return new ObitFArray, should be Unreffed when done
 */
ObitFArray* ObitConvUtilGaus (ObitImage *inImage, ofloat Beam[3])
{
  ObitFArray* outArray=NULL;
  ofloat      amp, Cen[2], GauMod[3];
  olong       ndim=2, naxis[2];

  /* Use Gaussian - make array from Input */
  naxis[0] = inImage->myDesc->inaxes[0];  
  naxis[1] = inImage->myDesc->inaxes[1]; 
  outArray = ObitFArrayCreate("Convolving function", ndim, naxis);
  
  /* Zero array */
  ObitFArrayFill (outArray, 0.0);
  
  /* Insert Gaussian */
  Cen[0] = (ofloat)(naxis[0]/2);      /* zero ref */
  Cen[1] = (ofloat)(naxis[1]/2);
  if (fabs(inImage->myDesc->cdelt[0])>0.0) {
    GauMod[0] = Beam[0] / fabs(inImage->myDesc->cdelt[0]); 
    GauMod[1] = Beam[1] / fabs(inImage->myDesc->cdelt[0]);
  } else { /* Beam in pixels */
    GauMod[0] = Beam[0];
    GauMod[1] = Beam[1];
  }
  GauMod[2] = Beam[2] - 90.0;
  /* Want Gaussian with normalized area */
  amp = 1.0 / (2.0 * G_PI * (GauMod[0]/2.3548) * (GauMod[1]/2.3548));
  
  /* Make elliptical Gaussian in outArray */
  ObitFArray2DEGauss (outArray, amp, Cen, GauMod);

  return outArray;
} /* end ObitConvUtilGaus */

/**
 * Deconvolves a Gaussian "beam" from a Gaussian component.  
 * Routine translated from the AIPSish APL/SUB/DECONV.FOR/DECONV  
 * \param fmaj    Convolved major axis 
 * \param fmin    Convolved minor axis 
 * \param fpa     Convolved position angle of major axis 
 * \param cmaj    Beam major axis 
 * \param cmin    Beam minor axis 
 * \param cpa     Beam position angle of major axis 
 * \param rmaj    [out] Actual major axis; = 0 => unable to fit 
 * \param rmin    [out] Actual  minor axis; = 0 => unable to fit 
 * \param rpa     [out] Actual position angle of major axis 
 */
void ObitConvUtilDeconv (ofloat fmaj, ofloat fmin, ofloat fpa, 
			 ofloat cmaj, ofloat cmin, ofloat cpa, 
			 ofloat *rmaj, ofloat *rmin, ofloat *rpa)
{
  ofloat      cmj2, cmn2, fmj2, fmn2, sinc, cosc, rhoc, 
    sigic2, det, rhoa, lfpa, lcpa, konst = 28.647888;
  olong csux;

  /* Get useful constants */
  csux = (olong) ((fpa+900.0)/180.0);
  lfpa = (fpa+900.0) - csux*180.0;
  csux = (olong) ((cpa+900.0)/180.0);
  lcpa = (cpa+900.0) - csux*180.0;

  cmj2 = cmaj * cmaj;
  cmn2 = cmin * cmin;
  fmj2 = fmaj * fmaj;
  fmn2 = fmin * fmin;
  sinc = (lfpa - lcpa) / konst;
  cosc = cos (sinc);
  sinc = sin (sinc);

  /* Trigonometry now */
  rhoc = (fmj2 - fmn2) * cosc - (cmj2 - cmn2);
  if (rhoc == 0.0) {
    sigic2 = 0.0;
    rhoa = 0.0;
  } else {
    sigic2 = atan((fmj2 - fmn2) * sinc / rhoc);
    rhoa = ((cmj2 - cmn2) - (fmj2 - fmn2) * cosc) / (2.0 * cos (sigic2));
  } 
  (*rpa) = sigic2 * konst + lcpa;
  det = ((fmj2 + fmn2) -(cmj2 + cmn2)) / 2.0;
  (*rmaj) = det - rhoa;
  (*rmin) = det + rhoa;

  /* Swap to get major > minor */
  (*rmaj) = MAX (0.0, *rmaj);
  (*rmin) = MAX (0.0, *rmin);
  (*rmaj) = sqrt (fabs (*rmaj));
  (*rmin) = sqrt (fabs (*rmin));
  if (*rmaj < *rmin) {
    sinc = (*rmaj);
    (*rmaj) = (*rmin);
    (*rmin) = sinc;
    (*rpa) = (*rpa)+90.0;
  } 

  /* Fix up PA */
  csux = (olong) ((*rpa+900.0)/180.0);
  *rpa = (*rpa+900.0) - csux*180.0;
  if (*rmaj == 0.0) {
    (*rpa) = 0.0;
  } else if (*rmin == 0.0) {
    if ((fabs(*rpa-lfpa) > 45.0)  &&  (fabs(*rpa-lfpa) < 135.0)) {
      csux = (olong) ((*rpa+450.0)/180.0);
      *rpa = (*rpa+450.0) - csux*180.0;
    }
  } 
} /* end of routine ObitConvUtilDeconv */ 
