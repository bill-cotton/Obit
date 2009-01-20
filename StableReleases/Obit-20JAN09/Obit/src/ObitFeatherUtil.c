/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2008                                          */
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
#include <string.h>
#include "ObitFeatherUtil.h"
#include "ObitImageUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFeatherUtil.c
 * ObitFeatherUtil utility function definitions for ObitImage class.
 *
 * Routines for feathering together images at different resolutions.
 */


/*----------------------Private functions---------------------------*/
/** Private: */

/*----------------------Public functions---------------------------*/
/**
 * Create an FFT for an efficient size equal to or larger than image
 * One needed for each direction to be FFTed.
 * \param in  Image to be FFTed
 * \param dir FFT direction (OBIT_FFT_Forward, OBIT_FFT_Reverse)
 * \return ObitFFT, must be Unrefed when done
 */
ObitFFT* ObitFeatherUtilCreateFFT (ObitImage *in, ObitFFTdir dir)
{
  ObitFFT *outFFT;
  olong effDim[2], rank = 2;
  gchar *name=NULL;

  /* Checks */
  g_assert(ObitImageIsA(in));

  /* Get image info from descriptor */
  effDim[0] = ObitFFTSuggestSize (in->myDesc->inaxes[0]);
  effDim[1] = ObitFFTSuggestSize (in->myDesc->inaxes[1]);

  name = g_strconcat ("FFT for ",in->name,NULL);
  outFFT    = newObitFFT(name, dir, OBIT_FFT_HalfComplex, rank, effDim);
  g_free(name);

  return outFFT;
} /* end ObitFeatherUtilCreateFFT */

/**
 * returns  half plane  CArray object of suitable size (2D) for FFTing image
 * \param inFFT   FFT to be applied
 * \return CArray, must be Unrefed when done
 */
ObitCArray* ObitFeatherUtilCreateFFTArray(ObitFFT *inFFT)
{
  ObitCArray *outCArray=NULL;
  olong ndim=2, naxis[2];

  /* Checks */
  g_assert(ObitFFTIsA(inFFT));

  /* FFT info */
  naxis[0] = 1 + inFFT->dim[0]/2;
  naxis[1] =  inFFT->dim[1];
  outCArray = ObitCArrayCreate ("Temp CArray for FFT", ndim, naxis);
  return outCArray;
} /* end  ObitFeatherUtilCreateFFTArray */

/**
 * Return scratch Image the size needed for a padded version of in
 * Returned image is 2-D
 * \param inFFT    FFT defining proper padded image size 
 *   (created by ObitFeatherUtilCreateFFT)
 * \param inImage  Obit Image to be padded.
 * \param err      ObitErr error stack.
 * \return Scratch Image, must be Unrefed when done
 */
ObitImage* ObitFeatherUtilCreatePadImage(ObitFFT *inFFT, 
					 ObitImage *inImage, ObitErr *err)
{
  ObitImage *out=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong size[IM_MAXDIM] = {1,1,1,1,1,1,1};
  gchar *routine = "ObitFeatherUtilCreatePadImage";

  /* Checks */
  g_assert(ObitFFTIsA(inFFT));
  g_assert(ObitImageIsA(inImage));

  /* Set dimension  - only need 2-D */
  dim[0] = 2;
  size[0] = (olong)inFFT->dim[0];
  size[1] = (olong)inFFT->dim[1];
  ObitInfoListAlwaysPut(inImage->info, "ScrSize", OBIT_long, dim, size);

  /* create Scratch */
  out = newObitImageScratch (inImage, err);
  if (err->error) Obit_traceback_val (err, routine, inImage->name, out);

  return out;
} /* end  ObitFeatherUtilCreatePadImage */

/**
 *  Zero Pads an image as needed for an FFT
 *  Any blanked values are replaced with zeroes
 * \param inFFT    FFT to be applied
 * \param inImage  Obit Image to be padded.
 * \param outImage Obit Image for output
 * \param err      ObitErr error stack.
 */
void ObitFeatherUtilPad (ObitFFT *inFFT, ObitImage *inImage, 
			 ObitImage *outImage, ObitErr *err)
{
  ObitFArray *outArray=NULL;
  olong ndim=2, naxis[2], pixOff[2];
  gchar *routine = "ObitFeatherUtilPad";

  /* Checks */
  g_assert(ObitFFTIsA(inFFT));
  g_assert(ObitImageIsA(inImage));
  g_assert(ObitImageIsA(outImage));
  if (err->error) return;

  /* Read input, Copy Descriptor */
  ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  ObitImageRead (inImage, NULL, err);
  ObitImageClose (inImage,err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
    
  /* FFT info */
  naxis[0] = inFFT->dim[0];
  naxis[1] = inFFT->dim[1];

  outArray = ObitFArrayCreate ("Padded array", ndim, naxis);

  /* Copy/pad/deblank array */
  ObitFArrayPad(inImage->image, outArray, 1.0);

  /* Copy selected descriptor info */
  outImage->myDesc = 
    ObitImageDescCopy (inImage->myDesc, outImage->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Reset output image size of first two axes */
  /* Update reference pixel, pixel shift an integral number */
  pixOff[0] = naxis[0]/2-inImage->myDesc->inaxes[0]/2;
  pixOff[1] = naxis[1]/2-inImage->myDesc->inaxes[1]/2;
  outImage->myDesc->crpix[0] = inImage->myDesc->crpix[0] + pixOff[0];
  outImage->myDesc->crpix[1] = inImage->myDesc->crpix[1] + pixOff[1];

  /* Update size */
  outImage->myDesc->naxis     = 2;
  outImage->myDesc->inaxes[0] = naxis[0];
  outImage->myDesc->inaxes[1] = naxis[1];
  outImage->myDesc->bitpix = -32; /* output floating */

  /* Write output image */
  ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  ObitImageWrite (outImage, outArray->array, err);
  ObitImageClose (outImage,err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Release pixel memory storage */
  inImage->image  = ObitFArrayUnref(inImage->image);
  outImage->image = ObitFArrayUnref(outImage->image);
  outArray = ObitFArrayUnref(outArray); /* release array */
} /* end ObitFeatherUtilPad */

/**
 * Increases the size of an image and zero pads, blanks replaced by zero
 * \param naxis    dimension array of desired output
 * \param inImage  Obit Image to be padded.
 * \param outImage Obit Image for output, Must previously exist
 * \param err      ObitErr error stack.
 */
void ObitFeatherUtilBigger (olong *naxis, ObitImage *inImage, 
			    ObitImage *outImage, ObitErr *err)
{
  ObitFArray *outArray=NULL;
  olong ndim=2, pixOff[2];
  gchar *routine = "ObitFeatherUtilBigger";

  /* Checks */
  g_assert(ObitImageIsA(inImage));
  g_assert(ObitImageIsA(outImage));
  g_assert(naxis!=NULL);
  if (err->error) return;

  /* Read input, Copy Descriptor */
  ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  ObitImageRead (inImage, NULL, err);
  ObitImageClose (inImage,err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
    
  outArray = ObitFArrayCreate ("Padded array", ndim, naxis);

  /* Reset output image size of first two axes */
  /* Update reference pixel, pixel shift an integral number */
  pixOff[0] = naxis[0]/2-inImage->myDesc->inaxes[0]/2;
  pixOff[1] = naxis[1]/2-inImage->myDesc->inaxes[1]/2;
  outImage->myDesc->crpix[0] = inImage->myDesc->crpix[0] + pixOff[0];
  outImage->myDesc->crpix[1] = inImage->myDesc->crpix[1] + pixOff[1];

  /* Update size */
  outImage->myDesc->inaxes[0] = naxis[0];
  outImage->myDesc->inaxes[1] = naxis[1];
  outImage->myDesc->bitpix = -32; /* output floating */

  /* Write output image */
  ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  ObitImageWrite (outImage, outArray->array, err);
  ObitImageClose (outImage,err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Release pixel memory storage */
  inImage->image  = ObitFArrayUnref(inImage->image);
  outImage->image = ObitFArrayUnref(outImage->image);
  outArray = ObitFArrayUnref(outArray); /* release array */
} /* end ObitFeatherUtilBigger */

/**
 * Zero fill array, any blanked values are replaced with zeroes
 * \param inFFT     Gives size of FFT needed
 * \param inArray   FArray to be padded.
 * \param outArray  FArray containing inArray but zero filled.
 */
void ObitFeatherUtilPadArray (ObitFFT *inFFT, ObitFArray *inArray, 
			      ObitFArray *outArray)
{
  olong pos1[2], pos2[2];

  /* Checks */
  g_assert(ObitFFTIsA(inFFT));
  g_assert(ObitFArrayIsA(inArray));
  g_assert(ObitFArrayIsA(outArray));

  /* Zero fill output */
  ObitFArrayFill(outArray, 0.0);

  /* Insert inArray into outArray - center as well as possible */
  pos1[0] = inFFT->dim[0]/2;     pos1[1] = inFFT->dim[1]/2;
  pos2[0] = inArray->naxis[0]/2; pos2[1] = inArray->naxis[1]/2;
  ObitFArrayShiftAdd (outArray, pos1, inArray, pos2, 1.0, outArray);

  /* Replace any blanks with zeroes */
  ObitFArrayDeblank(outArray, 0.0);
} /* end ObitFeatherUtilPadArray */

/**
 * Extract a Real array from one padded for FFTs
 * \param inFFT     Gives size of FFT needed
 * \param inArray   FArray with FFT results.
 * \param err   Obit error stack object.
 */
 ObitFArray* ObitFeatherUtilExtract (ObitFFT *inFFT, ObitFArray *inArray, 
				     ObitErr *err)
{
  ObitFArray *outArray=NULL;
  olong blc[2], trc[2], cen[2];
  gchar *routine = "ObitFeatherUtilExtract";

  /* Checks */
  g_assert(ObitFFTIsA(inFFT));
  g_assert(ObitFArrayIsA(inArray));
  if (err->error) return outArray;

  /* Get window to extract */
  cen[0] = inFFT->dim[0]/2;  cen[1] = inFFT->dim[1]/2; 
  blc[0] = cen[0] - inArray->naxis[0] / 2; trc[0] = cen[0] - 1 + inArray->naxis[0] / 2;
  blc[1] = cen[1] - inArray->naxis[1] / 2; trc[1] = cen[1] - 1 + inArray->naxis[1] / 2;

  /* Extract */
  outArray = ObitFArraySubArr(inArray, blc, trc, err);
  if (err->error) Obit_traceback_val (err, routine, inArray->name, outArray);

  return outArray;
} /* end ObitFeatherUtilExtract */

/**
 * Make uv plane weighting array
 * Creates an FArray the size of a plane in inImage, FFT,
 * takes real part and normalizes the central value to one
 * Resulting array is returned.
 * \param inImage  Obit Image to be padded.
 * \param inFFT    Obit forward FFT object
 * \param err      ObitErr error stack.
 * \return uv mask array
 */
ObitFArray* ObitFeatherUtilMakeBeamMask (ObitImage *inImage, ObitFFT *inFFT, 
					 ObitErr *err)
{
  olong ndim=2, naxis[2], pos[2];
  ofloat peak, norm;
  ObitFArray *maskArray=NULL, *outArray=NULL, *FFTArray=NULL;
  ObitCArray *uvArray=NULL;
  gchar *name=NULL;
  gchar *routine = "ObitFeatherUtilMakeBeamMask";

  /* Checks */
  g_assert(ObitImageIsA(inImage));
  g_assert(ObitFFTIsA(inFFT));
  if (err->error) return maskArray;

  /* Make copy of data array */
  naxis[0] = inImage->myDesc->inaxes[0]; 
  naxis[1] = inImage->myDesc->inaxes[1]; 
  outArray = ObitFArrayCreate("Beam mask", ndim, naxis);

  /* Read input */
  ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  ObitImageRead (inImage, NULL, err);
  ObitImageClose (inImage,err);
  if (err->error) Obit_traceback_val (err, routine, inImage->name, outArray);

  /* Add model */
  ObitFeatherUtilCreateModel(inImage, outArray);
  inImage->image = ObitFArrayUnref(inImage->image);  /* free memory */
    
  /* Pad for FFT */
  naxis[0] = inFFT->dim[0];  naxis[1] = inFFT->dim[1]; 
  FFTArray =  ObitFArrayCreate ("FFT array", ndim, naxis);
  ObitFeatherUtilPadArray (inFFT, outArray, FFTArray);
  outArray = ObitFArrayUnref(outArray);  /* Cleanup */

  /* Swaparoonie */
  ObitFArray2DCenter (FFTArray);

  /* FFT */
  uvArray = ObitFeatherUtilCreateFFTArray(inFFT);
  ObitFFTR2C (inFFT, FFTArray, uvArray);
  FFTArray = ObitFArrayUnref(FFTArray);  /* Cleanup */
  
  /* Extract Real part */
  naxis[0] = uvArray->naxis[0]; naxis[1] = uvArray->naxis[1];
  name = g_strconcat ("Mask array for ", inImage->name, NULL);
  maskArray = ObitFArrayCreate (name, ndim, naxis);
  g_free(name);

  ObitCArrayReal(uvArray, maskArray);
  uvArray = ObitCArrayUnref(uvArray); /* Cleanup */

  /* Normalize */
  pos[0] = 0; pos[1] = 1+naxis[1]/2;
  peak = ObitFArrayMax(maskArray, pos);
  if (peak!=0.0) norm = 1.0 / peak;
  else norm = 1.0;
  ObitFArraySMul(maskArray, norm);

  return maskArray;
} /* end ObitFeatherUtilMakeBeamMask */

/**
 * Fill an FArray with a model the size and shape of the resolution in an image.
 * A model image is inserted in outArray derived from the restoring beam in
 * image.  The size and geometry of outArray must be those described by image
 * \param image     Image with info, image member MUST be present
 * \param outArray  Python FArray to receive model, must previously exist.
 */
void ObitFeatherUtilCreateModel (ObitImage *image, ObitFArray *outArray)
{
  ofloat beamMaj, beamMin, beamPA, *cdelt, *crpix, *crota, amp, Cen[2], GauMod[3];
  olong *inaxes;

  /* Checks */
  g_assert(ObitImageIsA(image));
  g_assert(ObitFArrayIsA(outArray));
  g_assert(ObitFArrayIsCompatable(image->image, outArray));

  /* Get image info from descriptor */
  beamMaj = image->myDesc->beamMaj;
  beamMin = image->myDesc->beamMin;
  beamPA  = image->myDesc->beamPA;
  cdelt   = image->myDesc->cdelt;
  crpix   = image->myDesc->crpix;
  crota   = image->myDesc->crota;
  inaxes  = image->myDesc->inaxes;

  /* Check that beam OK */
  g_assert(beamMaj> 0.0001/3600.0);

  /* Zero array */
  ObitFArrayFill (outArray, 0.0);

  amp = 1.0;
  Cen[0] = (ofloat)(inaxes[0]/2);      /* zero ref */
  Cen[1] = (ofloat)(inaxes[1]/2);
  GauMod[0] = beamMaj/fabs(cdelt[0]);  /* Gaussian */
  GauMod[1] = beamMin/fabs(cdelt[0]);
  GauMod[2] = beamPA-90.0;

  /* Make elliptical Gaussian in ouArray */
  ObitFArray2DEGauss (outArray, amp, Cen, GauMod);
} /* end ObitFeatherUtilCreateModel */

/**
 * Accumulate the weighted FT of an FArray
 * inImage is FFTed, multiplied by wtArray and accumulated into accArray
 * \param FFTfor    FFT object to FT input array (on inImage)
 * \param inImage   Image to be accumulated must be a size compatable with 
 *                  FFTfor, returned with contents swapped for FFTs
 * \param wtArray   FArray containing accumulation weights, must be a size 
 *                  compatable with FT of input array.
 * \param accArray  CArray in which the results are to be accumulated.
 *                  must be a size compatable with FT of input array.
 * \param workArray CArray for temporary storage of FT of input array.
 * \param err     ObitErr error stack.
 */
void ObitFeatherUtilAccumImage (ObitFFT *FFTfor, ObitImage *inImage, 
				ObitFArray *wtArray, ObitCArray *accArray, 
				ObitCArray *workArray, ObitErr *err)
{
  ofloat beamMaj, beamMin, *cdelt, factor;
  gchar *routine = "ObitFeatherUtilAccumImage";

  /* Checks */
  g_assert(ObitFFTIsA(FFTfor));
  g_assert(ObitImageIsA(inImage));
  g_assert(ObitFArrayIsA(wtArray));
  g_assert(ObitCArrayIsA(accArray));
  g_assert(ObitCArrayIsA(workArray));
  if (err->error) return;

  /* Check compatability */
  Obit_return_if_fail(ObitCArrayIsCompatable (accArray, workArray), err,
		      "workArray, accArray incompatible");

  /* Read input */
  ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  ObitImageRead (inImage, NULL, err);
  ObitImageClose (inImage,err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Swaparoonie */
  ObitFArray2DCenter (inImage->image);
    
  /* FFT */
  ObitFFTR2C (FFTfor, inImage->image, workArray);

  inImage->image = ObitFArrayUnref(inImage->image); /* release array */

  /* Multiply by weights */
  ObitCArrayFMul(workArray, wtArray, workArray);

  /* Scale by inverse of beam area to get units the same */
  /* Get image info from descriptor */
  /* Python version multiplies beam? by 3600.0 */
  beamMaj = inImage->myDesc->beamMaj;
  beamMin = inImage->myDesc->beamMin;
  cdelt   = &inImage->myDesc->cdelt[0];
  factor = (fabs(cdelt[1])/beamMaj) * (fabs(cdelt[1])/beamMin);
  ObitCArraySMul(workArray, factor);

  /* Accumulate */
  ObitCArrayAdd(accArray, workArray, accArray);

} /* end ObitFeatherUtilAccumImage */

/**
 *  HGEOM-like operation (Before EWG got to it)
 *  Extract a Real array from one padded for FFTs
 * \param inImage   Image to be interpolated
 * \param tmplImage Image whose geometry is to be used.
 * \param outImage  values from inImage on grid of tmplImage
 *                  undefined values set to zero to allow FFT.
 * \param err       ObitErr error stack.
 */
void ObitFeatherUtilInterpol (ObitImage *inImage, ObitImage *tmplImage, 
			      ObitImage *outImage, ObitErr *err)
{
  olong i;
  olong inPlane[5] = {1,1,1,1,1};   /* Might need these passed */
  olong outPlane[5] = {1,1,1,1,1};
  olong hwidth = 2;
  gchar *routine = "ObitFeatherUtilInterpol";

  /* Checks */
  g_assert(ObitImageIsA(inImage));
  g_assert(ObitImageIsA(tmplImage));
  g_assert(ObitImageIsA(outImage));
  if (err->error) return;

  /* Clone, interpolate, deblank */
  ObitImageClone(tmplImage, outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Open output */
  ObitImageOpen (outImage, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Copy selected header info */
  outImage->myDesc->beamMaj = inImage->myDesc->beamMaj;
  outImage->myDesc->beamMin = inImage->myDesc->beamMin;
  outImage->myDesc->beamPA  = inImage->myDesc->beamPA;
  for (i=0; i<IMLEN_VALUE; i++) {
    outImage->myDesc->obsdat[i] = inImage->myDesc->obsdat[i];
    outImage->myDesc->teles[i]  = inImage->myDesc->teles[i];
    outImage->myDesc->instrument[i]  = inImage->myDesc->instrument[i];
    outImage->myDesc->observer[i]    = inImage->myDesc->observer[i];
  } 
  
  /* Force header update */
  outImage->myStatus = OBIT_Modified;

  /* Close output */
  ObitImageClose (outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  
  /* Interpolate image */
  ObitImageUtilInterpolateImage(inImage, outImage, inPlane, 
				outPlane, hwidth, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  
  /* reread image, deblank, rewrite */
  ObitImageOpen (outImage, OBIT_IO_ReadWrite, err);
  ObitImageRead (outImage, NULL, err);
  ObitDataIOSet ((ObitData*)outImage,err);
  ObitFArrayDeblank (outImage->image, 0.0);
  ObitImageWrite (outImage, NULL, err);
  ObitImageClose (outImage,err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Free pixel memory */
  outImage->image  = ObitFArrayUnref(outImage->image);
}
  /* end ObitFeatherUtilInterpol */

/**
 * Extract the subimage in inArray corresponding to outImage
 * This assumes that both input and output images have pixels
 * on the same locations (i.e. one the padded version of the other)
 * outImage is updated in permanent storage (disk)
 * \param inImage  Image describing inArray
 * \param inArray  array from which values are to be extracted
 * \param outImage accepts values from inImage, must exist and be fully defined
 * \param err     ObitErr error stack.
 */
void ObitFeatherUtilSubImage (ObitImage *inImage, ObitFArray *inArray, 
			      ObitImage *outImage, ObitErr *err)
{
  ObitFArray *outArray=NULL;
  olong blc[MAXFARRAYDIM] = {1,1,1,1,1,1,1};
  olong trc[MAXFARRAYDIM] = {0,0,0,0,0,0,0};
  gchar *routine = "ObitFeatherUtilSubImag";

  /* Checks */
  g_assert(ObitImageIsA(inImage));
  g_assert(ObitFArrayIsA(inArray));
  g_assert(ObitImageIsA(outImage));
  if (err->error) return;

  /* determine window in image */
  blc[0] = ((olong)(inImage->myDesc->crpix[0]) - (olong)(outImage->myDesc->crpix[0]));
  blc[1] = ((olong)(inImage->myDesc->crpix[1]) - (olong)(outImage->myDesc->crpix[1]));
  
  trc[0] = blc[0] + outImage->myDesc->inaxes[0] - 1;
  trc[1] = blc[1] + outImage->myDesc->inaxes[1] - 1;
  
  /* Extract array */
  outArray = ObitFArraySubArr(inArray, blc, trc, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  ObitImageWrite (outImage, outArray->array, err);
  ObitImageClose (outImage,err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Release pixel memory storage */
  outImage->image = ObitFArrayUnref(outImage->image);
  outArray = ObitFArrayUnref(outArray); /* release array */
} /* end  ObitFeatherUtilSubImage */
