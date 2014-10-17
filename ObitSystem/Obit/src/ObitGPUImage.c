/* $Id:  $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014                                               */
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
#include "ObitGPUImage.h"
#include "ObitGPUFArray.h"
#include "ObitGPUFInterpolate.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGPUImage.c
 *
 * ObitGPUImage GPU enhanced Image utilities
 *
 * Image utililties using a GPU.
 * 
 */
/*---------------Public functions---------------------------*/
/**
 * Fill the pixels in outImage by interpolation to the corresponding 
 * locations in inImage.
 * There is no interpolation between planes
 * GPU Implementation
 * The GPU should be assigned and reset externally to this routine.
 * \param finterp  GPU Interpolator
 * \param inImage  Image to be interpolated
 * \param outImage Image to be written.  Must be previously instantiated.
 * \param inPlane  desired plane in inImage, 1-rel pixel numbers on planes 3-7
 * \param outPlane desired plane in outImage
 * \param hwidth   interpolation halfwidth (1 or 2 usually OK, 4 max)
 * \param err      Error stack, returns if not empty.
 */
void 
ObitGPUImageInterpolateImageXY (ObitGPUFInterpolate *finterp, 
				ObitImage *inImage, ObitImage *outImage, 
				olong *inPlane, olong *outPlane,
				ObitErr *err)
{
  ObitIOSize IOBy;
  ObitImageDesc *tmpDesc=NULL;
  olong iblc[IM_MAXDIM], itrc[IM_MAXDIM], oblc[IM_MAXDIM], otrc[IM_MAXDIM];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, j;
  odouble RAPnt, DecPnt;
  ofloat fblank = ObitMagicF();
  gchar *today=NULL;
  gchar *routine = "ObitGPUImageInterpolateImageXY";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));
  g_assert (inPlane!=NULL);
  g_assert (outPlane!=NULL);

  for (i=0; i<IM_MAXDIM; i++) iblc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) itrc[i] = 0;
  for (i=0; i<IM_MAXDIM; i++) oblc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) otrc[i] = 0;

  /* Do I/O by plane and all of plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListPut (inImage->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
  ObitInfoListPut (outImage->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
  /* Get any previous blc, trc */
  ObitInfoListGetTest (inImage->info, "BLC", &type, dim, iblc); 
  ObitInfoListGetTest (inImage->info, "TRC", &type, dim, itrc);
  dim[0] = 7;
  for (i=0; i<5; i++) iblc[i+2] = itrc[i+2] = inPlane[i];
  ObitInfoListPut (inImage->info, "BLC", OBIT_long, dim, iblc, err); 
  ObitInfoListPut (inImage->info, "TRC", OBIT_long, dim, itrc, err);
  for (i=0; i<5; i++) oblc[i+2] = otrc[i+2] = outPlane[i];
  ObitInfoListPut (outImage->info, "BLC", OBIT_long, dim, oblc, err); 
  ObitInfoListPut (outImage->info, "TRC", OBIT_long, dim, otrc, err);

  /* Open images */
  if ((ObitImageOpen (inImage, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, inImage->name);
    return;
  }
  if ((ObitImageOpen (outImage, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, outImage->name);
    return;
  }
  /* Adjust output descriptor on first plane - copy from input */
  if ((outPlane[0]==1) && (outPlane[1]==1) && (outPlane[2]==1) && (outPlane[3]==1) 
      && (outPlane[4]==1)) {
    /* Copy of old descriptor */
    tmpDesc = ObitImageDescCopy (outImage->myDesc, tmpDesc, err);
    /* update Descriptive stuff from input */
    ObitImageDescCopyDesc (inImage->myDesc, outImage->myDesc, err);
    if (err->error) Obit_traceback_msg (err, routine, inImage->name);

    /* Creation date today */
    today = ObitToday();
    strncpy (outImage->myDesc->date, today, IMLEN_VALUE-1);
    if (today) g_free(today);

    /* Precess pointing position if necessary */
    if (inImage->myDesc->equinox!=tmpDesc->equinox) {
      ObitImageDescGetPoint (inImage->myDesc, &RAPnt, &DecPnt);
      if ((fabs(inImage->myDesc->equinox-1950.0)<0.01) && 
	  (fabs(tmpDesc->equinox-2000.0)<0.01))
	ObitSkyGeomBtoJ (&RAPnt, &DecPnt);
      else if ((fabs(inImage->myDesc->equinox-2000.0)<0.01) && 
	       (fabs(tmpDesc->equinox-1950.0)<0.01))
	ObitSkyGeomJtoB (&RAPnt, &DecPnt);
      outImage->myDesc->obsra  = RAPnt;
      outImage->myDesc->obsdec = DecPnt;
    }

    /* restore first two planes geometry */
    outImage->myDesc->epoch   = tmpDesc->epoch;
    outImage->myDesc->equinox = tmpDesc->equinox;
    for (j=0; j<2; j++) {
      outImage->myDesc->inaxes[j] = tmpDesc->inaxes[j];
      outImage->myDesc->cdelt[j]  = tmpDesc->cdelt[j];
      outImage->myDesc->crota[j]  = tmpDesc->crota[j];
      outImage->myDesc->crpix[j]  = tmpDesc->crpix[j];
      outImage->myDesc->crval[j]  = tmpDesc->crval[j];
      for (i=0; i<IMLEN_KEYWORD; i++) outImage->myDesc->ctype[j][i] = tmpDesc->ctype[j][i];
    }
    tmpDesc = ObitImageDescUnref(tmpDesc);
  }
  
  /* Read input plane */
  if ((ObitImageRead (inImage,NULL , err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR reading image %s", 
		   routine, inImage->name);
    return;
  }
  /* Convert pure zero to fblank */
  ObitFArrayInClip (inImage->image, -1.0e-25, 1.0e-25, fblank);

  /* Interpolate */
  ObitGPUFInterpolateImage(finterp, inImage->image, outImage->image, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Write output */
  if ((ObitImageWrite (outImage, NULL, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR writing image %s", 
		   routine, outImage->name);
    return;
  }

  /* Close */
  if ((ObitImageClose (outImage, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, outImage->name);
    return;
  }
  if ((ObitImageClose (inImage, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, inImage->name);
    return;
  }

  if ((ObitImageClose (outImage, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, outImage->name);
    return;
  }
} /* end  ObitGPUImageInterpolateImageXY */

