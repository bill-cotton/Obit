/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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

#include "ObitPolnUnwind.h"
#include "ObitThread.h"
#include "ObitSinCos.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif /* VELIGHT */
/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/*----------------------Private functions---------------------------*/
static void Unwrap (ObitFArray* RMFArray, odouble Lamb2, 
		    ObitFArray* inQFArray,  ObitFArray* inUFArray, 
		    ObitFArray* outQFArray, ObitFArray* outUFArray);
/*----------------------Public functions---------------------------*/
/**
 * Unrotate the Faraday rotation in a set of Q, U cubes.
 * If the input images are ObitImageMF like images then the first plane 
 * of the output images will be the average of the corrected images.
 * \param rmImage   Image with rotation measures to be removed in the 
 *                  first plane
 * \param inQImage  Q Image cube to be unwound
 * \param inUImage  U Image cube to be unwound
 * \param outQImage output derotated Q image
 * \param outUImage output derotated U image
 * \param err       Obit error stack object.
 */
void ObitPolnUnwindCube (ObitImage *rmImage, ObitImage *inQImage, ObitImage *inUImage, 
			 ObitImage *outQImage, ObitImage *outUImage, ObitErr *err)
{
  olong nplane, iplane, plane[5]={1,1,1,1,1}, noffset=0;
  olong naxis[2], nx, ny;
  ofloat norm=1.0;
  ObitIOSize IOBy;
  ObitInfoType type;
  ObitIOCode retCode;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean accum=FALSE;
  odouble freq, Lamb2;
  ObitFArray *RMFArray=NULL, *inQFArray=NULL, *inUFArray=NULL, 
    *outQFArray=NULL, *outUFArray=NULL, 
    *sumQFArray=NULL, *sumUFArray=NULL;
  gchar *today=NULL, keyword[9];
  gchar *routine = "ObitPolnUnwindCube";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail((ObitImageIsA(inQImage)), err, "%s: inQImage not an Image", 
			routine);
  Obit_return_if_fail((ObitImageIsA(inUImage)), err, "%s: inUImage not an Image", 
			routine);
  Obit_return_if_fail((ObitImageIsA(rmImage)), err, "%s: rmImage not an Image", 
			routine);

  /* Open input images to get info */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListAlwaysPut (inQImage->info, "IOBy", OBIT_long, dim, &IOBy);
  inQImage->extBuffer = TRUE;   /* Using inFArrays as I/O buffer */
  retCode = ObitImageOpen (inQImage, OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, inQImage->name);
  ObitInfoListAlwaysPut (inUImage->info, "IOBy", OBIT_long, dim, &IOBy);
  inUImage->extBuffer = TRUE;   /* Using inFArrays as I/O buffer */
  retCode = ObitImageOpen (inUImage, OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, inUImage->name);
  ObitInfoListAlwaysPut (rmImage->info, "IOBy", OBIT_long, dim, &IOBy);
  rmImage->extBuffer = TRUE;   /* Using inFArrays as I/O buffer */
  retCode = ObitImageOpen (rmImage, OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, rmImage->name);

  /* Check compatability */
  Obit_return_if_fail(((inQImage->myDesc->inaxes[0]==inUImage->myDesc->inaxes[0]) && 
		       (inQImage->myDesc->inaxes[1]==inUImage->myDesc->inaxes[1]) &&
		       (inQImage->myDesc->inaxes[2]==inUImage->myDesc->inaxes[2])), err,
		      "%s: Input images incompatable", routine);
  Obit_return_if_fail(((rmImage->myDesc->inaxes[0]==inUImage->myDesc->inaxes[0]) && 
		       (rmImage->myDesc->inaxes[1]==inUImage->myDesc->inaxes[1])), err,
		      "%s: Input/RM images incompatable", routine);

  /* What plane does spectral data start on */
  if (!strncmp(inQImage->myDesc->ctype[inQImage->myDesc->jlocf], "SPECLNMF", 8)) {
    noffset = 2;  /* What plane does spectral data start on */
    ObitInfoListGetTest (inQImage->myDesc->info, "NTERM", &type, dim, &noffset);
    accum = TRUE;   /* Accumulate sum */
  } else {   /* Normal spectral cube */
    noffset = 0;
    accum = FALSE;   /* Do not accumulate sum */
  }

  /* Determine number of frequency planes */
  nplane = inQImage->myDesc->inaxes[inQImage->myDesc->jlocf]-noffset;

 /* Create image buffers */
  nx = inQImage->myDesc->inaxes[0];
  ny = inQImage->myDesc->inaxes[1];
  naxis[0] = nx;  naxis[1] = ny; 
  RMFArray   = ObitFArrayCreate (NULL, 2, naxis);
  inQFArray  = ObitFArrayCreate (NULL, 2, naxis);
  inUFArray  = ObitFArrayCreate (NULL, 2, naxis);
  outQFArray = ObitFArrayCreate (NULL, 2, naxis);
  outUFArray = ObitFArrayCreate (NULL, 2, naxis);
  if (accum) {  /* Accumulate sum? */
   sumQFArray = ObitFArrayCreate (NULL, 2, naxis);
   sumUFArray = ObitFArrayCreate (NULL, 2, naxis);
   /* Zero Fill */
   ObitFArrayFill (sumQFArray, 0.0);
   ObitFArrayFill (sumUFArray, 0.0);
 }

  /* Clone output from input */
  ObitImageClone (inQImage, outQImage, err);
  ObitImageClone (inUImage, outUImage, err);
  if (err->error) goto cleanup;

  /* Output Image descriptors */
  outQImage->myDesc = ObitImageDescCopy (inQImage->myDesc, outQImage->myDesc, err);
  outUImage->myDesc = ObitImageDescCopy (inUImage->myDesc, outUImage->myDesc, err);
  if (err->error) goto cleanup;

  /* Creation date today */
  today = ObitToday();
  strncpy (outQImage->myDesc->date, today, IMLEN_VALUE);
  strncpy (outUImage->myDesc->date, today, IMLEN_VALUE);
  if (today) g_free(today);

  /* Open Output images */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListAlwaysPut (outQImage->info, "IOBy", OBIT_long, dim, &IOBy);
  outQImage->extBuffer = TRUE;   /* Using outFArrays as I/O buffer */
  retCode = ObitImageOpen (outQImage, OBIT_IO_ReadWrite, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;
  ObitInfoListAlwaysPut (outUImage->info, "IOBy", OBIT_long, dim, &IOBy);
  outUImage->extBuffer = TRUE;   /* Using outFArrays as I/O buffer */
  retCode = ObitImageOpen (outUImage, OBIT_IO_ReadWrite, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Zero fill output */
  for (iplane=0; iplane<outQImage->myDesc->inaxes[2]; iplane++) {
    plane[0] = iplane+1;
    retCode = ObitImagePutPlane (outQImage, outQFArray->array, plane, err);
    retCode = ObitImagePutPlane (outUImage, outUFArray->array, plane, err);
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;
  } /* end loop over planes */

  /* Read RM plane */
  plane[0] = 1;  /* Select correct plane */
  retCode = ObitImageGetPlane (rmImage, RMFArray->array, plane, err);
  /* if it didn't work bail out */
  if ((retCode!=OBIT_IO_OK) || (err->error))goto cleanup; 
  retCode = ObitImageClose (rmImage, err);
  rmImage->extBuffer = FALSE;   /* May need I/O buffer later */

  /* Loop over Q, U, planes */
  for (iplane=0; iplane<nplane; iplane++) {
    /* Lambda^2 Check for MFImage outputs */
    if (!strncmp(inQImage->myDesc->ctype[inQImage->myDesc->jlocf], "SPECLNMF", 8)) {
	sprintf (keyword, "FREQ%4.4d",iplane+1);
	freq = inQImage->myDesc->crval[inQImage->myDesc->jlocf] + 
	  inQImage->myDesc->cdelt[inQImage->myDesc->jlocf] * 
	  (inQImage->myDesc->plane - inQImage->myDesc->crpix[inQImage->myDesc->jlocf]);
	ObitInfoListGetTest (inQImage->myDesc->info, keyword, &type, dim, &freq);
       } else {   /* Normal spectral cube */
	freq = inQImage->myDesc->crval[inQImage->myDesc->jlocf] + 
	  inQImage->myDesc->cdelt[inQImage->myDesc->jlocf] * 
	  (inQImage->myDesc->plane - inQImage->myDesc->crpix[inQImage->myDesc->jlocf]);
      }
    Lamb2 = (VELIGHT/freq)*(VELIGHT/freq);

    plane[0] = iplane+noffset+1;  /* Select correct plane */
    retCode = ObitImageGetPlane (inQImage, inQFArray->array, plane, err);
    retCode = ObitImageGetPlane (inUImage, inUFArray->array, plane, err);
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

    /* Correct to output */
    Unwrap(RMFArray, Lamb2, inQFArray, inUFArray, outQFArray, outUFArray);

    if (accum) {  /* Accumulate sum? */
      ObitFArrayAdd (sumQFArray, outQFArray, sumQFArray);
      ObitFArrayAdd (sumUFArray, outUFArray, sumUFArray);
    }

    /* Write output planes */
    plane[0] = iplane+noffset+1;  /* Select correct plane */
    retCode = ObitImagePutPlane (outQImage, outQFArray->array, plane, err);
    retCode = ObitImagePutPlane (outUImage, outUFArray->array, plane, err);
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;
  } /* end loop over planes */

  if (accum) {  /* Write accumulated Q, U images to first plane */
    /* Normalize */
    norm = 1.0 / nplane;
    ObitFArraySMul (sumQFArray, norm);
    ObitFArraySMul (sumUFArray, norm);
    plane[0] = 1;
    retCode = ObitImagePutPlane (outQImage, sumQFArray->array, plane, err);
    retCode = ObitImagePutPlane (outUImage, sumUFArray->array, plane, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) goto cleanup;
  }

  /* Cleanup */
 cleanup:
  /* Close Q,U inputs */
  retCode = ObitImageClose (inQImage, err);
  inQImage->extBuffer = FALSE;   /* May need I/O buffer later */
  retCode = ObitImageClose (inUImage, err);
  inUImage->extBuffer = FALSE;   /* May need I/O buffer later */

  /* Close Q,U outputs */
  retCode = ObitImageClose (outQImage, err);
  outQImage->extBuffer = FALSE;   /* May need I/O buffer later */
  retCode = ObitImageClose (outUImage, err);
  outUImage->extBuffer = FALSE;   /* May need I/O buffer later */

  /* Delete arrays */
  RMFArray   = ObitFArrayUnref(RMFArray);
  inQFArray  = ObitFArrayUnref(inQFArray); 
  inUFArray  = ObitFArrayUnref(inUFArray);
  outQFArray = ObitFArrayUnref(outQFArray);
  outUFArray = ObitFArrayUnref(outUFArray);
  sumQFArray = ObitFArrayUnref(sumQFArray);
  sumUFArray = ObitFArrayUnref(sumUFArray);
  if (err->error) Obit_traceback_msg (err, routine, outQImage->name);
} /* end ObitPolnUnwindCube */

/**
 * Remove rotation measure from Q, U, pixels 
 * \param RMFArray    Array of RMs (rad/m**2)
 * \param Lamb2       Lambda**2 in m**2
 * \param inQFArray   Input Q pixels
 * \param inUFArray   Input U pixels
 * \param outQFArray  Output Q pixels
 * \param outUFArray  Output U pixels
 */
static void Unwrap (ObitFArray* RMFArray, odouble Lamb2, 
		    ObitFArray* inQFArray, ObitFArray* inUFArray, 
		    ObitFArray* outQFArray, ObitFArray* outUFArray)
{
  olong i;
  ofloat q, u, s, c, arg, fblank = ObitMagicF();

  for (i=0; i<RMFArray->arraySize; i++) {
    q = inQFArray->array[i];
    u = inUFArray->array[i];
    if ((RMFArray->array[i]!=fblank) && (q!=fblank) && 	(u!=fblank)) {
      arg = -2.0*RMFArray->array[i] * Lamb2;
      ObitSinCosCalc(arg, &s, &c);
      outQFArray->array[i] = q * c - u * s;
      outUFArray->array[i] = q * s + u * c;
   } else {   /* Blanked */
      outQFArray->array[i] = fblank;
      outUFArray->array[i] = fblank;
    }
  }  /* End loop over array */
} /* end  Unwrap */

