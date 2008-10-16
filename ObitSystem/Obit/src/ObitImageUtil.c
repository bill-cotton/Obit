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

#include <sys/types.h>
#include <time.h>
#include "ObitImageUtil.h"
#include "ObitUVWeight.h"
#include "ObitSkyGeom.h"
#include "ObitFInterpolate.h"
#include "ObitPBUtil.h"
#include "ObitIOImageFITS.h"
#include "ObitFArray.h"
#include "ObitFArrayUtil.h"
#include "ObitConvUtil.h"
#include "ObitFeatherUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageUtil.c
 * ObitImageUtil module function definitions for ObitImage class.
 */

/*---------------Private function prototypes----------------*/
/** Private: Get Date string for current date */
static void ObitImageUtilCurDate (gchar *date, olong len);

/** Private: Fit Beam size to dirty beam */
static void ObitImageUtilFitBeam (ObitImage *beam, ObitErr *err);

/*----------------------Public functions---------------------------*/

/**
 * Create basic ObitImage structure and fill out descriptor.
 * Imaging parameters are on the inUV info member as arrays for a number 
 * of fields.
 * \li "nChAvg" OBIT_int (1,1,1) number of channels to average.
 *              This is for spectral line observations and is ignored
 *              if the IF axis on the uv data has more than one IF.
 *              Default is continuum = average all freq/IFs. 0=> all.
 * \li "rotate" OBIT_float (?,1,1) Desired rotation on sky (from N thru E) in deg. [0]
 * \li "nx"     OBIT_int (?,1,1) Dimension of image in RA [no default].
 *              This and the following are arrays with one entry per field.
 * \li "nxBeam" OBIT_int (?,1,1) Dimension of beam in RA, [def. nx]
 * \li "ny"     OBIT_int (?,1,1) Dimension of image in declination[no default]
 * \li "nyBeam" OBIT_int (?,1,1) Dimension of beam in declination, [def. ny]
 * \li "xCells" OBIT_float (?,1,1) X (=RA) cell spacing in asec [no default]
 * \li "yCells" OBIT_float (?,1,1) Y (=dec) cell spacing in asec [no default]
 * \li "xShift" OBIT_float (?,1,1) Desired shift in X (=RA) in degrees. [0]
 * \li "yShift" OBIT_float (?,1,1) Desired shift in Y (=dec) in degrees. [0]
 * \param inUV     Input uv data. 
 * \param fieldNo  Which field (1-rel) in imaging parameter arrays.
 * \param doBeam   if TRUE also create beam as the myBeam member of 
 *                 returned image.
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created ObitImage.
 */
ObitImage* ObitImageUtilCreateImage (ObitUV *inUV, olong fieldNo, 
				     gboolean doBeam, ObitErr *err)
{
  ObitImage *outImage=NULL, *theBeam=NULL;
  gchar outName[121];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong *iarray, nChAvg, nx=0, ny=0, nxBeam=0, nyBeam=0;
  ofloat *farray, xCells=0.0, yCells=0.0, rotate, xShift, yShift;
  gchar *routine = "ObitImageUtilCreateImage";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outImage;
  g_assert (ObitUVIsA(inUV));
  g_assert (fieldNo>0);

  /* open/close uv data to fully instantiate if not already open */
  if (inUV->myStatus==OBIT_Inactive) {
    ObitUVFullInstantiate (inUV, TRUE, err);
    if (err->error) Obit_traceback_val (err, routine, inUV->name, outImage);
  }

  /* frequency tables - always refresh in case selection redefined */
  ObitUVGetFreq (inUV, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outImage);

 /* Create output - name to include field name */
  g_snprintf (outName, 120, "%s Field  %d",inUV->name, fieldNo);
  outImage = newObitImage(outName);

  /* Need beam as well? */
  if (doBeam) {
    g_snprintf (outName, 120, "%s Beam  %d",inUV->name, fieldNo);
    theBeam = newObitImage(outName);

    /* save the beam on output */
    outImage->myBeam = (Obit*)theBeam;
  } /* end create beam */

  /* Get parameters for image */
  /* Number of channels to average, defaults to all */
  nChAvg = 1000000000; /* unlikely number to exceed */
  ObitInfoListGetTest(inUV->info, "nChAvg", &type, dim, (gpointer)&nChAvg);
  if (nChAvg<=0) nChAvg = 1000000000;

  /* Image size */
  if (!ObitInfoListGetP(inUV->info, "nx", &type, dim, (gpointer*)&iarray)) {
    Obit_log_error(err, OBIT_Error, 
		   "ObitImageUtilCreateImage: %s MUST define nx", 
		   inUV->name);
  }
  /* Make sure fields large enough */
  if ((!err->error) && (dim[0]<fieldNo)) {
    Obit_log_error(err, OBIT_Error, 
		   "ObitImageUtilCreateImage: %s Parameter array too small", 
		   inUV->name);
    return outImage;
  }
  if(iarray!=NULL) nx = iarray[fieldNo-1];

  if (!ObitInfoListGetP(inUV->info, "ny", &type, dim, (gpointer*)&iarray)) {
    Obit_log_error(err, OBIT_Error, 
		   "ObitImageUtilCreateImage: %s MUST define ny", 
		   inUV->name);
  }
  if(iarray!=NULL) ny = iarray[fieldNo-1];

  /* Beam size */
  if (doBeam) {
    if (!ObitInfoListGetP(inUV->info, "nxBeam", &type, dim, (gpointer*)&iarray)) {
      nxBeam = nx; /* not found - use default */
    } else {
      if (iarray!=NULL) nxBeam = iarray[fieldNo-1]; /* passed value */
    }
    if (!ObitInfoListGetP(inUV->info, "nyBeam", &type, dim, (gpointer*)&iarray)) {
      nyBeam = ny; /* not found - use default */
    } else {
      if (iarray!=NULL) nyBeam = iarray[fieldNo-1]; /* passed value */
    }
  } /* end beam size */

  /* Cell Spacing */
  if (!ObitInfoListGetP(inUV->info, "xCells", &type, dim, (gpointer*)&farray)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: %s MUST define xCells", 
		   routine, inUV->name);
  }
  if (farray!=NULL) xCells = farray[fieldNo-1]/3600.0;

  if (!ObitInfoListGetP(inUV->info, "yCells", &type, dim, (gpointer*)&farray)) {
    Obit_log_error(err, OBIT_Error, 
		   "ObitImageUtilCreateImage: %s MUST define yCells", 
		   inUV->name);
  }
  if (farray!=NULL) yCells = farray[fieldNo-1]/3600.0;

  /* rotation  default 0 */
  rotate = 0.0;
  ObitInfoListGetP(inUV->info, "rotate", &type, dim, (gpointer*)&farray);
  if (farray!=NULL) rotate = farray[fieldNo-1]; /* passed value */

  /* Field shift */
  xShift = 0.0;
  yShift = 0.0;
  ObitInfoListGetP(inUV->info, "xShift", &type, dim, (gpointer*)&farray);
  if (farray!=NULL) xShift = farray[fieldNo-1];

  ObitInfoListGetP(inUV->info, "yShift", &type, dim, (gpointer*)&farray);
  if (farray!=NULL) yShift = farray[fieldNo-1];

  /* bail out if an error so far */
  if (err->error) return outImage;

  /* Set values on descriptor(s) */
  outImage->myDesc->xshift = xShift;
  outImage->myDesc->yshift = yShift;
  outImage->myDesc->crota[0] = 0.0;
  outImage->myDesc->crota[1] = rotate;
  outImage->myDesc->cdelt[0] = xCells;
  outImage->myDesc->cdelt[1] = yCells;
  outImage->myDesc->inaxes[0] = nx;
  outImage->myDesc->inaxes[1] = ny;

  /* Temporarily add image to shift on uv data */
  xShift = inUV->myDesc->xshift;
  yShift = inUV->myDesc->yshift;
  inUV->myDesc->xshift += outImage->myDesc->xshift;
  inUV->myDesc->yshift += outImage->myDesc->yshift;

  /* Fill in descriptor */
  ObitImageUtilUV2ImageDesc (inUV->myDesc, outImage->myDesc, nChAvg);

  /* Header may have changed */
  if (outImage->myStatus!=OBIT_Inactive) outImage->myStatus = OBIT_Modified;

  /* Also beam if needed */
  if (doBeam) {
   theBeam->myDesc->xshift = xShift;
   theBeam->myDesc->yshift = yShift;
   theBeam->myDesc->crota[0] = 0.0;
   theBeam->myDesc->crota[1] = rotate;
   theBeam->myDesc->cdelt[0] = xCells;
   theBeam->myDesc->cdelt[1] = yCells;
   theBeam->myDesc->inaxes[0] = nxBeam;
   theBeam->myDesc->inaxes[1] = nyBeam;
   /* Fill in descriptor */
   ObitImageUtilUV2ImageDesc (inUV->myDesc, theBeam->myDesc, nChAvg);
   /* Header may have changed */
   if (theBeam->myStatus!=OBIT_Inactive) theBeam->myStatus = OBIT_Modified;
 }

  /* Replace uv data shift */
  inUV->myDesc->xshift = xShift;
  inUV->myDesc->yshift = yShift;

  return outImage;
} /* end ObitImageUtilCreateImage */

/**
 * Grids, FFTs and makes corrections for the gridding convolution.
 * Uses (creating if necessary) the myGrid member of out.
 * \param inUV     Input uv data. Should be in form of stokes to be imaged
 *                 will all calibration and selection applied and any 
 *                 weighting applied.
 * \param outImage Image to be written.  Must be previously instantiated.
 *                 Beam normalization factor is written to output Beam
 *                 infoList as SUMWTS
 * \param doBeam   if TRUE also make beam.  Will make the myBeam member of 
 *                 outImage.
 *                 If FALSE, and myGrid->BeamNorm 0.0 then reads SUMWTS value 
 *                 from beam infolist
 * \param doWeight if TRUE Apply uniform weighting corrections to uvdata
 *                 before imaging.
 * Weighting parameters on inUV:
 * \li "nuGrid" OBIT_long scalar = Number of "U" pixels in weighting grid.
 *              [defaults to "nx"]
 * \li "nvGrid" OBIT_long scalar = Number of "V" pixels in weighting grid.
 * \li "WtBox"  OBIT_long scalar = Size of weighting box in cells [def 1]
 * \li "WtFunc" OBIT_long scalar = Weighting convolution function [def. 1]
 *              1=Pill box, 2=linear, 3=exponential, 4=Gaussian
 *              if positive, function is of radius, negative in u and v.
 * \li "xCells" OBIT_float scalar = Image cell spacing in X in asec.
 * \li "yCells" OBIT_float scalar = Image cell spacing in Y in asec.
 * \li "UVTaper" OBIT_float scalar = UV taper width in kilowavelengths. [def. no taper].
 *              NB: If the taper is applied her is should not also be applied
 *              in the imaging step as the taper will be applied to the
 *              output data.
 * \li "Robust" OBIT_float scalar = Briggs robust parameter. [def. 0.0]
 *              < -7 -> Pure Uniform weight, >7 -> Pure natural weight.
 *              Uses AIPS rather than Briggs definition of Robust.
 * \li "WtPower" OBIT_float scalar = Power to raise weights to.  [def = 1.0]
 *              Note: a power of 0.0 sets all the output weights to 1 as modified
 *              by uniform/Tapering weighting.  Applied in determinng weights 
 *              as well as after.
 * \param channel  Which frequency channel to image, 0->all.
 * \param err      Error stack, returns if not empty.
 */
void ObitImageUtilMakeImage (ObitUV *inUV, ObitImage *outImage, 
			     olong channel, gboolean doBeam, gboolean doWeight,
			     ObitErr *err)
{
  ObitImage *theBeam=NULL;
  ObitIOSize IOBy;
  ObitInfoType type;
  ofloat sumwts, imMax, imMin;
  gchar *outName=NULL;
  olong i, ichannel, icLo, icHi, pos[5];
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "ObitImageUtilMakeImage";

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitImageIsA(outImage));
  if (channel > outImage->myDesc->inaxes[outImage->myDesc->jlocf]) {
    Obit_log_error(err, OBIT_Error, 
		   "Requested channel %d > max %d for %s", 
		   channel, 
		   outImage->myDesc->inaxes[outImage->myDesc->jlocf], 
		   outImage->name);
    return;
 }

  /* Apply uniform weighting? */
  if (doWeight) ObitUVWeightData (inUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* Need new gridding member? */
  outName = g_strconcat ("UVGrid for: ",inUV->name,NULL);
  if (outImage->myGrid == NULL) outImage->myGrid = newObitUVGrid(outName);
  g_free(outName);

  /*  Open image ReadOnly to get proper descriptor */
  dim[0] = 7;
  for (i=0; i<IM_MAXDIM; i++) {blc[i] = 1; trc[i] = 0;}
  ObitInfoListPut (outImage->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (outImage->info, "TRC", OBIT_long, dim, trc, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  if ((ObitImageOpen (outImage, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, 
		   "ERROR opening image %s", outImage->name);
    return;
  }
  ObitImageClose (outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* Get beam normalization factor from beam if needed */
  if (!doBeam) {
    /* Get Beam member from outImage */
    theBeam = (ObitImage*)outImage->myBeam;
    g_assert (ObitImageIsA(theBeam));

    ObitInfoListGet(theBeam->info, "SUMWTS", &type, dim, (gpointer)&sumwts, err);
    if (err->error) Obit_traceback_msg (err, routine, theBeam->name);
    outImage->myGrid->BeamNorm = sumwts;
  }    

  /* Loop over channels selected */
  icLo = 1; 
  icHi = outImage->myDesc->inaxes[outImage->myDesc->jlocf];
  if (channel>0) {icLo = channel; icHi = channel;}
  for (ichannel=icLo; ichannel<=icHi; ichannel++) {

    /* Make beam if needed */
    if (doBeam) {
      /* Get Beam member from outImage */
      theBeam = (ObitImage*)outImage->myBeam;
      g_assert (ObitImageIsA(theBeam));
      
      /* Set blc, trc */
      blc[2] = ichannel;
      trc[0] = outImage->myDesc->inaxes[0];
      trc[1] = outImage->myDesc->inaxes[1];
      trc[2] = ichannel;
      
      IOBy = OBIT_IO_byPlane;
      dim[0] = 1;
      ObitInfoListPut (theBeam->info, "IOBy", OBIT_long, dim,
		       (gpointer)&IOBy, err);
      dim[0] = 7;
      ObitInfoListPut (theBeam->info, "BLC", OBIT_long, dim,
		       (gpointer)blc, err); 
      ObitInfoListPut (theBeam->info, "TRC", OBIT_long, dim,
		       (gpointer)trc, err);

      /*  Open beam to get descriptor */
      if ((ObitImageOpen (theBeam, OBIT_IO_ReadWrite, err) 
	   != OBIT_IO_OK) || (err->error>0)) { /* error test */
	Obit_log_error(err, OBIT_Error, 
		       "ERROR opening image %s", theBeam->name);
      }
      
      /* Gridding setup */
      ObitUVGridSetup (outImage->myGrid, inUV, theBeam->myDesc, 
		       outImage->myDesc, TRUE, err);
      if (err->error) Obit_traceback_msg (err, routine, theBeam->name);
      
      /* Grid Beam */
      ObitUVGridReadUV (outImage->myGrid, inUV, err);
      if (err->error) Obit_traceback_msg (err, routine, theBeam->name);
      
      /* FFT, Gridding correction */
      ObitUVGridFFT2Im(outImage->myGrid, theBeam->image, err);
      if (err->error) Obit_traceback_msg (err, routine, theBeam->name);
      
      /* Write image */
      ObitImageWrite (theBeam, NULL, err);
      if (err->error) Obit_traceback_msg (err, routine, theBeam->name);
      
      /* Tell Sum of gridding weights */
      Obit_log_error(err, OBIT_InfoErr, 
		     "Sum of Weights %g for %s",outImage->myGrid->BeamNorm,outImage->name);

      /* tell Max/Min */
      imMax = ObitFArrayMax (theBeam->image, pos);
      imMin = ObitFArrayMin (theBeam->image, pos);
      Obit_log_error(err, OBIT_InfoErr, 
		     "Beam max %g, min %g for %s", imMax, imMin, theBeam->name);

      /* Fit beam */
      ObitImageUtilFitBeam (theBeam, err);
      if (err->error) Obit_traceback_msg (err, routine, theBeam->name);

      /* Close Image */
      ObitImageClose (theBeam, err);
      if (err->error) Obit_traceback_msg (err, routine, theBeam->name);

      /* Free image buffer */
      theBeam->image = ObitFArrayUnref(theBeam->image);

   } /* end making beam */
    
    /* Now make image */
    /* Set blc, trc */
    blc[2] = ichannel;
    trc[0] = outImage->myDesc->inaxes[0];
    trc[1] = outImage->myDesc->inaxes[0];
    trc[2] = ichannel;

    IOBy = OBIT_IO_byPlane;
    dim[0] = 1;
    ObitInfoListPut (outImage->info, "IOBy", OBIT_long, dim,
		     (gpointer)&IOBy, err);
    dim[0] = 7;
    ObitInfoListPut (outImage->info, "BLC", OBIT_long, dim,
		     (gpointer)blc, err); 
    ObitInfoListPut (outImage->info, "TRC", OBIT_long, dim,
		     (gpointer)trc, err);

    /*  Open image to get descriptor */
    if ((ObitImageOpen (outImage, OBIT_IO_ReadWrite, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, 
		     "ERROR opening image %s", outImage->name);
    }
    

    /* Gridding setup */
    ObitUVGridSetup (outImage->myGrid, inUV, theBeam->myDesc, 
		     outImage->myDesc, FALSE, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    
    /* Grid Image */
    ObitUVGridReadUV (outImage->myGrid, inUV, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    
    /* Gridding correction */
    ObitUVGridFFT2Im(outImage->myGrid, outImage->image, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    
    /* Write image */
    ObitImageWrite (outImage, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    
    /* tell Max/Min */
    imMax = ObitFArrayMax (outImage->image, pos);
    imMin = ObitFArrayMin (outImage->image, pos);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Image max %g, min %g for %s", imMax, imMin,  outImage->name);
    /* ObitErrLog(err);  Progress Report */
    
    /* Copy Gaussian beam fit */
    outImage->myDesc->beamMaj = theBeam->myDesc->beamMaj;
    outImage->myDesc->beamMin = theBeam->myDesc->beamMin;
    outImage->myDesc->beamPA  = theBeam->myDesc->beamPA;
    
    /* Save last beam normalization in Beam infoList as "SUMWTS" */
    sumwts = outImage->myGrid->BeamNorm;
    dim[0] = 1;
    dim[1] = 0;
    ObitInfoListPut(theBeam->info, "SUMWTS", OBIT_float, dim, (gpointer)&sumwts, err);
    if (err->error) Obit_traceback_msg (err, routine, theBeam->name);

    /* Make sure Stokes correct */
    outImage->myDesc->crval[outImage->myDesc->jlocs] = 
      inUV->myDesc->crval[inUV->myDesc->jlocs];
    
    /* Close Image */
    ObitImageClose (outImage, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    
    /* Free image buffer */
    outImage->image = ObitFArrayUnref(outImage->image);

    /* Free gridding member */
    outImage->myGrid = ObitUVGridUnref(outImage->myGrid);
    
  } /* end loop over channels */
}  /* end ObitImageUtilMakeImage */

/**
 * Fill the pixels in outImage by interpolation to the corresponding locations
 * in inImage.
 * There is no interpolation between planes
 * \param inImage  Image to be interpolated.
 * \param outImage Image to be written.  Must be previously instantiated.
 * \param inPlane  desired plane in inImage, 1-rel pixel numbers on planes 3-7
 * \param outPlane desired plane in outImage
 * \param hwidth   interpolation halfwidth (1 or 2 usually OK, 4 max)
 * \param err      Error stack, returns if not empty.
 */
void 
ObitImageUtilInterpolateImage (ObitImage *inImage, ObitImage *outImage, 
			       olong *inPlane, olong *outPlane,
			       olong hwidth, ObitErr *err)
{
  ObitIOSize IOBy;
  ObitFInterpolate *interp=NULL;
  ObitImageDesc *tmpDesc=NULL;
  olong iblc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong itrc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong oblc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong otrc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, j, ix, iy, indx, pos[2];
  ofloat inPixel[2], outPixel[2], *out, fblank =  ObitMagicF();
  odouble RAPnt, DecPnt;
  gboolean OK;
  gchar *today=NULL;
  gchar *routine = "ObitImageUtilInterpolateImage";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));
  g_assert (inPlane!=NULL);
  g_assert (outPlane!=NULL);
 
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

  /* Make interpolator */
  interp = newObitFInterpolateCreate ("Interpolator", inImage->image, 
				      inImage->myDesc, hwidth);

  /* Get output aray pointer */
  pos[0] = pos[1] = 0;
  out = ObitFArrayIndex (outImage->image, pos);

  /* Loop over image interpolating */
  for (iy = 1; iy<=outImage->myDesc->inaxes[1]; iy++) { /* loop in y */
    outPixel[1] = (ofloat)iy;
    for (ix = 1; ix<=outImage->myDesc->inaxes[0]; ix++) {/* loop in x */
      outPixel[0] = (ofloat)ix;

      /* Get pixel in input image */
      OK = ObitImageDescCvtPixel (outImage->myDesc, inImage->myDesc, outPixel, inPixel, err);
      if (err->error) {
	interp = ObitFInterpolateUnref(interp); /* Free interpolator */
	Obit_traceback_msg (err, routine, inImage->name);
      }

      indx = (iy-1) * outImage->myDesc->inaxes[0] + (ix-1);
      if (OK) { /* In image? */
	/* interpolate */
	out[indx] = ObitFInterpolatePixel (interp, inPixel, err);
	if (err->error) Obit_traceback_msg (err, routine, inImage->name);
      } else { /* not in input */
	out[indx] = fblank;
      }

    } /* end loop over x */
  } /* end loop over y */
  

  /* Free interpolator */
  interp = ObitFInterpolateUnref(interp);

  /* Write output */
  if ((ObitImageWrite (outImage, NULL, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR writing image %s", 
		   routine, outImage->name);
    return;
  }

  /* Close */
  if ((ObitImageClose (inImage, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, inImage->name);
    return;
  }
  /* Free image buffer if not memory resident */
  if (inImage->mySel->FileType!=OBIT_IO_MEM) 
    inImage->image = ObitFArrayUnref(inImage->image);

  if ((ObitImageClose (outImage, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, outImage->name);
    return;
  }
  /* Free image buffer if not memory resident */
  if (outImage->mySel->FileType!=OBIT_IO_MEM) 
    outImage->image = ObitFArrayUnref(outImage->image);
} /* end  ObitImageUtilInterpolateImage */

/**
 * Fill the pixels in outImage by interpolation to the corresponding locations
 * in inImage given a Zernike model of distortions in in.
 * There is no interpolation between planes
 * \param inImage  Image to be interpolated.
 * \param outImage Image to be written.  Must be previously instantiated.
 * \param inPlane  desired plane in inImage, 1-rel pixel numbers on planes 3-7
 * \param outPlane desired plane in outImage
 * \param hwidth   interpolation halfwidth (1 or 2 usually OK, 4 max)
 * \param nZern    Number of Zernike terms, can handle up to 17
 * \param ZCoef    Array of Zernike coefficients (piston ignored)
 * \param err      Error stack, returns if not empty.
 */
void 
ObitImageUtilInterpolateImageZern (ObitImage *inImage, ObitImage *outImage, 
				   olong *inPlane, olong *outPlane,
				   olong hwidth, olong nZern, ofloat *ZCoef, 
				   ObitErr *err)
{
  ObitIOSize IOBy;
  ObitFInterpolate *interp=NULL;
  ObitImageDesc *tmpDesc=NULL;
  olong iblc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong itrc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong oblc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong otrc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, j, ix, iy, indx, pos[2];
  ofloat inPixel[2], outPixel[2], *out, fblank =  ObitMagicF();
  gboolean OK;
  gchar *today=NULL;
  gchar *routine = "ObitImageUtilInterpolateImageZern";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));
  g_assert (inPlane!=NULL);
  g_assert (outPlane!=NULL);
 
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
  ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  ObitImageOpen (outImage, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

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
  ObitImageRead (inImage,NULL , err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Make interpolator */
  interp = newObitFInterpolateCreate ("Interpolator", inImage->image, 
				      inImage->myDesc, hwidth);

  /* Get output aray pointer */
  pos[0] = pos[1] = 0;
  out = ObitFArrayIndex (outImage->image, pos);

  /* Loop over image interpolating */
  for (iy = 1; iy<=outImage->myDesc->inaxes[1]; iy++) { /* loop in y */
    outPixel[1] = (ofloat)iy;
    for (ix = 1; ix<=outImage->myDesc->inaxes[0]; ix++) {/* loop in x */
      outPixel[0] = (ofloat)ix;

      /* Get pixel in input image */
      OK = ObitImageDescCvtZern (outImage->myDesc, inImage->myDesc, nZern, ZCoef, 
				 outPixel, inPixel, err);
      if (err->error) Obit_traceback_msg (err, routine, inImage->name);

      indx = (iy-1) * outImage->myDesc->inaxes[0] + (ix-1);
      if (OK) { /* In image? */
	/* interpolate */
	out[indx] = ObitFInterpolatePixel (interp, inPixel, err);
	if (err->error) Obit_traceback_msg (err, routine, inImage->name);
      } else { /* not in input */
	out[indx] = fblank;
      }

    } /* end loop over x */
  } /* end loop over y */
  

  /* Free interpolator */
  interp = ObitFInterpolateUnref(interp);

  /* Write output */
  ObitImageWrite (outImage, NULL, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Close */
  ObitImageClose (inImage, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Free image buffer if not memory resident */
  if (inImage->mySel->FileType!=OBIT_IO_MEM) 
    inImage->image = ObitFArrayUnref(inImage->image);

  ObitImageClose (outImage, err) ;
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Free image buffer if not memory resident */
  if (outImage->mySel->FileType!=OBIT_IO_MEM) 
    outImage->image = ObitFArrayUnref(outImage->image);
} /* end  ObitImageUtilInterpolateImageZern */

/**
 * Fill the pixels in outImage by interpolation to the corresponding locations
 * in inImage.
 * Also calculates a weight based on a circle defined by radius; this is 1.0 
 * in the center and tapers with distance^2 to 0.0 outside.
 * If memOnly then the input image plane is assumed in inImage and only memory
 * resident parts of outImage and outWeight are modified.
 * There is no interpolation between planes
 * \param inImage   Image to be interpolated.
 * \param outImage  Image to be written.  Must be previously instantiated.
 * \param outWeight Weight image to be written.  Must be previously instantiated and
 *                  have same geometry as outImage.
 * \param memOnly   if TRUE then work only in memory
 * \param radius    Radius in pixels of weighting circle
 * \param inPlane   Desired plane in inImage, 1-rel pixel numbers on planes 3-7; 
 *                  ignored if memOnly
 * \param outPlane  Desired plane in outImage; ignored if memOnly
 * \param hwidth    Interpolation halfwidth (1 or 2 usually OK, 4 max)
 * \param err       Error stack, returns if not empty.
 */
void 
ObitImageUtilInterpolateWeight (ObitImage *inImage, ObitImage *outImage, 
				ObitImage *outWeight, gboolean memOnly,
				gint radius, olong *inPlane, olong *outPlane,
				olong hwidth, ObitErr *err)
{
  ObitIOSize IOBy;
  ObitFInterpolate *interp=NULL;
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong ix, iy, indx, pos[2];
  ofloat inPixel[2], outPixel[2], *out, *outWt, rad2, dist2, irad2;
  ofloat *crpix, wt, val, fblank =  ObitMagicF();
  gboolean OK;
  gchar *routine = "ObitImageUtilInterpolateImage";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));
  g_assert (ObitImageIsA(outWeight));
  g_assert (inPlane!=NULL);
  g_assert (outPlane!=NULL);
 
  /* Working version of radius */
  rad2 = radius * radius;
  irad2 = 1.0 / rad2;

  /* Do I/O by plane and all of plane */
  if (!memOnly) {
    IOBy = OBIT_IO_byPlane;
    dim[0] = 1;
    ObitInfoListPut (inImage->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
    ObitInfoListPut (outImage->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
    ObitInfoListPut (outWeight->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
    dim[0] = IM_MAXDIM;
    for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = trc[i+2] = inPlane[i];
    ObitInfoListPut (inImage->info, "BLC", OBIT_long, dim, blc, err); 
    ObitInfoListPut (inImage->info, "TRC", OBIT_long, dim, trc, err);
    for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = trc[i+2] = outPlane[i];
    ObitInfoListPut (outImage->info, "BLC", OBIT_long, dim, blc, err); 
    ObitInfoListPut (outImage->info, "TRC", OBIT_long, dim, trc, err);
    ObitInfoListPut (outWeight->info, "BLC", OBIT_long, dim, blc, err); 
    ObitInfoListPut (outWeight->info, "TRC", OBIT_long, dim, trc, err);

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
    
    if ((ObitImageOpen (outWeight, OBIT_IO_ReadWrite, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		     routine, outWeight->name);
      return;
    }
    /* Read input plane */
    if ((ObitImageRead (inImage,NULL , err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR reading image %s", 
		     routine, inImage->name);
      return;
    }
  } /* end of not memory only */

  /* Make interpolator */
  interp = newObitFInterpolateCreate ("Interpolator", inImage->image, 
				      inImage->myDesc, hwidth);

  /* Get output aray pointer */
  pos[0] = pos[1] = 0;
  out   = ObitFArrayIndex (outImage->image, pos);
  outWt = ObitFArrayIndex (outWeight->image, pos);

  /* Coordinate reference pixel of input */
  crpix = inImage->myDesc->crpix;

  /* Loop over image interpolating */
  for (iy = 1; iy<=outImage->myDesc->inaxes[1]; iy++) { /* loop in y */
    outPixel[1] = (ofloat)iy;
    for (ix = 1; ix<=outImage->myDesc->inaxes[0]; ix++) {/* loop in x */
      outPixel[0] = (ofloat)ix;

     /* Get pixel in input image */
      OK = ObitImageDescCvtPixel (outImage->myDesc, inImage->myDesc, outPixel, inPixel, err);
      if (err->error) {
	interp = ObitFInterpolateUnref(interp);
	Obit_traceback_msg (err, routine, inImage->name);
      }

      if (OK) { /* In image? */
	/* weight based on distance from center of inImage */
	dist2 = (crpix[0]-inPixel[0])*(crpix[0]-inPixel[0]) + 
	  (crpix[1]-inPixel[1])*(crpix[1]-inPixel[1]);
	/*dist2 = (crpix[0]-xyzi[0])**2 + (crpix[1]-xyzi[1])**2;*/
	if (dist2 <= rad2) {
	  wt = 1.0 - dist2 * irad2;
	  wt = MAX (0.001, wt);
	} else {
	  wt = fblank;
	}
      } else {
	wt = fblank;  /* don't bother */
      }
      
      /* interpolate */
      /* array index in out for this pixel */
      indx = (iy-1) * outImage->myDesc->inaxes[0] + (ix-1);
      if (wt != fblank ) {
	val = ObitFInterpolatePixel (interp, inPixel, err);
	if (val != fblank ) val *= wt;
	out[indx] = val;
	if (err->error) {
	  interp = ObitFInterpolateUnref(interp);
	  Obit_traceback_msg (err, routine, inImage->name);
	}
	/* debug
	   if (out[indx] >1.0e5) fprintf (stderr,"bad int %d  %d %g %g \n",
					  ix,iy,out[indx], wt); */
	outWt[indx] = wt;
      } else {
	out[indx]   = fblank;
	outWt[indx] = fblank;
      }
      
    } /* end loop over x */
  } /* end loop over y */
  
  /* Free interpolator */
  interp = ObitFInterpolateUnref(interp);

  if (!memOnly) {
    /* Write output */
    if ((ObitImageWrite (outImage, NULL, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR writing image %s", 
		     routine, outImage->name);
      return;
    }
    
    /* Close */
    /* Write output */
    if ((ObitImageWrite (outWeight, NULL, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR writing image %s", 
		     routine, outWeight->name);
      return;
    }
    
    if ((ObitImageClose (inImage, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		     routine, inImage->name);
      return;
    }
    /* Free image buffer  if not memory resident */
    if (inImage->mySel->FileType!=OBIT_IO_MEM) 
      inImage->image = ObitFArrayUnref(inImage->image);
    if ((ObitImageClose (outImage, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		     routine, outImage->name);
      return;
    }
    /* Free image buffer if not memory resident */
    if (outImage->mySel->FileType!=OBIT_IO_MEM) 
      outImage->image = ObitFArrayUnref(outImage->image);
    if ((ObitImageClose (outWeight, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		     routine, outImage->name);
      return;
    }
    /* Free image buffer if not memory resident */
    if (outWeight->mySel->FileType!=OBIT_IO_MEM) 
      outWeight->image = ObitFArrayUnref(outWeight->image);
    
  } /* end of not memory only */
} /* end  ObitImageUtilInterpolateImage */

/**
 * Make antenna primary beam correction to an image based on the pointing
 * position in another image.
 * For frequencies < 1 GHz uses the VLA polynomial gain curves,
 * for higher frequencies, it uses a jinc function based on the antenna size.
 * \param inImage  Image to be corrected
 * \param pntImage Image with pointing position
 * \param outImage Image to be written.  Must be previously instantiated.
 * \param inPlane   Desired plane in inImage, 1-rel pixel numbers on planes 3-7; 
 *                  ignored if memOnly
 * \param outPlane  Desired plane in outImage; ignored if memOnly
 * \param antSize  Antenna size, used to correct beam for freq>1 GHz, def. 25m.
 * \param err      Error stack, returns if not empty.
 */
void 
ObitImageUtilPBCorr (ObitImage *inImage, ObitImage *pntImage, ObitImage *outImage, 
		     olong *inPlane, olong *outPlane, ofloat antSize, ObitErr *err)
{
  ObitIOSize IOBy;
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong ix, iy, indx, pos[2];
  ofloat inPixel[2], *in, *out;
  odouble RAPnt, DecPnt, Freq, ra, dec, xx, yy, zz, dist;
  ofloat pbf, equinox, fblank = ObitMagicF();
  gboolean doJinc, bad;
  ObitImageDesc *inDesc;
  gchar *routine = "ObitImageUtilPBCorr";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(pntImage));
  g_assert (ObitImageIsA(outImage));

  if (antSize<0.01) antSize = 25.0; /* default antenna size */

  /* Get pointing position */
  pntImage->extBuffer = TRUE;  /* Don't need buffer */
  if ((ObitImageOpen (pntImage, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, pntImage->name);
    return;
  }
  /* Use "Observed" position if given */
  RAPnt   = pntImage->myDesc->obsra;
  DecPnt  = pntImage->myDesc->obsdec;
  equinox = pntImage->myDesc->equinox ;
  if ((abs(RAPnt)<1.0e-5) && (abs(DecPnt)<1.0e-5)) {
    /* if zeroes - use reference position */
    RAPnt  = pntImage->myDesc->crval[pntImage->myDesc->jlocr];
    DecPnt = pntImage->myDesc->crval[pntImage->myDesc->jlocd];
  }
  RAPnt  *= DG2RAD;
  DecPnt *= DG2RAD;
  if ((ObitImageClose (pntImage, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, pntImage->name);
    return;
  }
  pntImage->extBuffer = FALSE;  /* May need buffer later */

  /* Do I/O by plane and all of plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListPut (inImage->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
  ObitInfoListPut (outImage->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
  dim[0] = 7;
  for (i=0; i<5; i++) blc[i+2] = trc[i+2] = inPlane[i];
  ObitInfoListPut (inImage->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (inImage->info, "TRC", OBIT_long, dim, trc, err);
  for (i=0; i<5; i++) blc[i+2] = trc[i+2] = outPlane[i];
  ObitInfoListPut (outImage->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (outImage->info, "TRC", OBIT_long, dim, trc, err);

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
  
  /* Read input plane */
  if ((ObitImageRead (inImage, NULL, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR reading image %s", 
		   routine, inImage->name);
    return;
  }

  /* Check that input and output are compatible */
  if (!ObitFArrayIsCompatable(inImage->image, outImage->image)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Input (%s) and output (%s) images are incompatible", 
		   routine, inImage->name, outImage->name);
    return;
  }

  /* check coordinate types */
  if (inImage->myDesc->coordType != pntImage->myDesc->coordType) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Input (%s) and pointing (%s) images have different coordinate types ", 
		   routine, inImage->name, pntImage->name);
    return;
  }

  /* Precess pointing position is needed */
  if ((abs(equinox-1950.0)<0.01) && 
      (abs(inImage->myDesc->equinox-2000.0)<0.01))
    ObitSkyGeomBtoJ (&RAPnt, &DecPnt);
  else if ((abs(equinox-2000.0)<0.01) && 
      (abs(inImage->myDesc->equinox-1950.0)<0.01))
    ObitSkyGeomJtoB (&RAPnt, &DecPnt);

  /* Get output aray pointer */
  pos[0] = pos[1] = 0;
  in  = ObitFArrayIndex (inImage->image, pos);
  out = ObitFArrayIndex (outImage->image, pos);

  /* Set up - get frequency, default 1 GHz*/
  Freq = 1.0e9;
  if (inImage->myDesc->jlocf>=0) 
    Freq = inImage->myDesc->crval[inImage->myDesc->jlocf];
  /* which beam model to use */
  doJinc = (Freq >= 1.0e9);
  inDesc = inImage->myDesc; /* Input descriptor */

  /* Loop over image  */
  for (iy = 1; iy<=inDesc->inaxes[1]; iy++) { /* loop in y */
    inPixel[1] = (ofloat)iy;
    for (ix = 1; ix<=inDesc->inaxes[0]; ix++) {/* loop in x */
      inPixel[0] = (ofloat)ix;

      /* array index in in and out for this pixel */
      indx = (iy-1) * inDesc->inaxes[0] + (ix-1);

      /* Is this pixel valid? */
      if (in[indx] != fblank) {
	/* Convert input pixel to position */
	bad = 
	  ObitSkyGeomWorldPos(inPixel[0], inPixel[1],
			      inDesc->crval[inDesc->jlocr], inDesc->crval[inDesc->jlocd],
			      inDesc->crpix[inDesc->jlocr], inDesc->crpix[inDesc->jlocd],
			      inDesc->cdelt[inDesc->jlocr], inDesc->cdelt[inDesc->jlocd],
			      inDesc->crota[inDesc->jlocd], &inDesc->ctype[inDesc->jlocr][4],
			      &ra, &dec);
	if (bad!=0) {
	  Obit_log_error(err, OBIT_Error, 
			 "%s: Error %d determining location of pixel in %s", 
			 routine, bad, inImage->name);
	  return;
	}
	
	/* Separation from pointing center */
	xx = DG2RAD * ra;
	yy = DG2RAD * dec;
	zz = sin (yy) * sin (DecPnt) + cos (yy) * cos (DecPnt) * cos (xx-RAPnt);
	zz = MIN (zz, 1.000);
	dist = acos (zz) * RAD2DG;
	
	/* primary beam correction */
	if (doJinc) {
	  pbf = ObitPBUtilJinc (dist, Freq, antSize, 0.0);
	} else {
	  pbf = ObitPBUtilPoly (dist, Freq, 0.0);
	} 

	/* Make correction */
	out[indx] = in[indx] / pbf;
      } else out[indx] = fblank; /* Don't correct if invalid */

    } /* end loop over x */
  } /* end loop over y */
  

  /* Write output */
  if ((ObitImageWrite (outImage, NULL, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR writing image %s", 
		   routine, outImage->name);
    return;
  }

  /* Close */
  if ((ObitImageClose (inImage, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, inImage->name);
    return;
  }
  /* Free image buffer if not memory resident */
  if (inImage->mySel->FileType!=OBIT_IO_MEM) 
    inImage->image = ObitFArrayUnref(inImage->image);
  if ((ObitImageClose (outImage, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, outImage->name);
    return;
  }
  /* Free image buffer  if not memory resident*/
  if (outImage->mySel->FileType!=OBIT_IO_MEM) 
    outImage->image = ObitFArrayUnref(outImage->image);

} /* end  ObitImageUtilPBCorr */

/**
 * Multiply antenna primary beam pattern by an image based on the pointing
 * position in another image.
 * For frequencies < 1 GHz uses the VLA polynomial gain curves,
 * for higher frequencies, it uses a jinc function based on the antenna size.
 * \param inImage  Image to be corrected
 * \param pntImage Image with pointing position
 * \param outImage Image to be written.  Must be previously instantiated.
 * \param inPlane   Desired plane in inImage, 1-rel pixel numbers on planes 3-7; 
 *                  ignored if memOnly
 * \param outPlane  Desired plane in outImage; ignored if memOnly
 * \param antSize  Antenna size
 * \param err      Error stack, returns if not empty.
 */
void 
ObitImageUtilPBApply (ObitImage *inImage, ObitImage *pntImage, ObitImage *outImage, 
		     olong *inPlane, olong *outPlane, ofloat antSize, ObitErr *err)
{
  ObitIOSize IOBy;
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong ix, iy, indx, pos[2];
  ofloat inPixel[2], *in, *out;
  odouble RAPnt, DecPnt, Freq, ra, dec, xx, yy, zz, dist ;
  ofloat pbf, equinox, fblank = ObitMagicF();
  gboolean doJinc, bad;
  ObitImageDesc *inDesc;
  gchar *routine = "ObitImageUtilPBApply";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(pntImage));
  g_assert (ObitImageIsA(outImage));

  if (antSize<0.01) antSize = 25.0; /* default antenna size */

  /* Get pointing position */
  pntImage->extBuffer = TRUE;  /* Don't need buffer */
  if ((ObitImageOpen (pntImage, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, pntImage->name);
    return;
  }

  /* Use  "Observed" position if given */
  RAPnt   = pntImage->myDesc->obsra;
  DecPnt  = pntImage->myDesc->obsdec;
  equinox = pntImage->myDesc->equinox ;
  if ((abs(RAPnt)<1.0e-5) && (abs(DecPnt)<1.0e-5)) {
    /* if zeroes - use reference position */
    RAPnt  = pntImage->myDesc->crval[pntImage->myDesc->jlocr];
    DecPnt = pntImage->myDesc->crval[pntImage->myDesc->jlocd];
  }
  RAPnt  *= DG2RAD;
  DecPnt *= DG2RAD;
  if ((ObitImageClose (pntImage, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, pntImage->name);
    return;
  }
  pntImage->extBuffer = FALSE;  /* May need buffer later */

  /* Do I/O by plane and all of plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListPut (inImage->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
  ObitInfoListPut (outImage->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
  dim[0] = 7;
  for (i=0; i<5; i++) blc[i+2] = trc[i+2] = inPlane[i];
  ObitInfoListPut (inImage->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (inImage->info, "TRC", OBIT_long, dim, trc, err);
  for (i=0; i<5; i++) blc[i+2] = trc[i+2] = outPlane[i];
  ObitInfoListPut (outImage->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (outImage->info, "TRC", OBIT_long, dim, trc, err);

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
  
  /* Read input plane */
  if ((ObitImageRead (inImage, NULL , err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR reading image %s", 
		   routine, inImage->name);
    return;
  }

  /* Check that input and output are compatible */
  if (!ObitFArrayIsCompatable(inImage->image, outImage->image)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Input (%s) and output (%s) images are incompatible", 
		   routine, inImage->name, outImage->name);
    return;
  }

  /* check coordinate types */
  if (inImage->myDesc->coordType != pntImage->myDesc->coordType) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Input (%s) and pointing (%s) images have different coordinate types ", 
		   routine, inImage->name, pntImage->name);
    return;
  }

  /* Precess pointing position is needed */
  if ((abs(equinox-1950.0)<0.01) && 
      (abs(inImage->myDesc->equinox-2000.0)<0.01))
    ObitSkyGeomBtoJ (&RAPnt, &DecPnt);
  else if ((abs(equinox-2000.0)<0.01) && 
      (abs(inImage->myDesc->equinox-1950.0)<0.01))
    ObitSkyGeomJtoB (&RAPnt, &DecPnt);

  /* Get output aray pointer */
  pos[0] = pos[1] = 0;
  in  = ObitFArrayIndex (inImage->image, pos);
  out = ObitFArrayIndex (outImage->image, pos);

  /* Set up - get frequency, default 1 GHz*/
  Freq = 1.0e9;
  if (inImage->myDesc->jlocf>=0) 
    Freq = inImage->myDesc->crval[inImage->myDesc->jlocf];
  /* which beam model to use */
  doJinc = (Freq >= 1.0e9);
  inDesc = inImage->myDesc; /* Input descriptor */

  /* Loop over image  */
  for (iy = 1; iy<=inDesc->inaxes[1]; iy++) { /* loop in y */
    inPixel[1] = (ofloat)iy;
    for (ix = 1; ix<=inDesc->inaxes[0]; ix++) {/* loop in x */
      inPixel[0] = (ofloat)ix;

      /* array index in in and out for this pixel */
      indx = (iy-1) * inDesc->inaxes[0] + (ix-1);

      /* Is this pixel valid? */
      if (in[indx] != fblank) {
	/* Convert input pixel to position */
	bad = 
	  ObitSkyGeomWorldPos(inPixel[0], inPixel[1],
			      inDesc->crval[inDesc->jlocr], inDesc->crval[inDesc->jlocd],
			      inDesc->crpix[inDesc->jlocr], inDesc->crpix[inDesc->jlocd],
			      inDesc->cdelt[inDesc->jlocr], inDesc->cdelt[inDesc->jlocd],
			      inDesc->crota[inDesc->jlocd], &inDesc->ctype[inDesc->jlocr][4],
			      &ra, &dec);
	if (bad!=0) {
	  Obit_log_error(err, OBIT_Error, 
			 "%s: Error %d determining location of pixel in %s", 
			 routine, bad, inImage->name);
	  return;
	}
	
	/* Separation from pointing center */
	xx = DG2RAD * ra;
	yy = DG2RAD * dec;
	zz = sin (yy) * sin (DecPnt) + cos (yy) * cos (DecPnt) * cos (xx-RAPnt);
	zz = MIN (zz, 1.000);
	dist = acos (zz) * RAD2DG;
	
	/* primary beam correction */
	if (doJinc) {
	  pbf = ObitPBUtilJinc (dist, Freq, antSize, 0.0);
	} else {
	  pbf = ObitPBUtilPoly (dist, Freq, 0.0);
	} 

	/* debug
	if (pbf>0.9) {
	  fprintf (stderr,"pdf %f @ %d  %d %lf\n",pbf,ix,iy,dist* RAD2DG*3600.0);
	} */
 
	/* Make correction */
	out[indx] = in[indx] * pbf;
      } else out[indx] = fblank; /* Don't correct if invalid */

    } /* end loop over x */
  } /* end loop over y */
  

  /* Write output */
  if ((ObitImageWrite (outImage, NULL, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR writing image %s", 
		   routine, outImage->name);
    return;
  }

  /* Close */
  if ((ObitImageClose (inImage, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, inImage->name);
    return;
  }
  /* Free image buffer if not memory resident */
  if (inImage->mySel->FileType!=OBIT_IO_MEM) 
    inImage->image = ObitFArrayUnref(inImage->image);

  if ((ObitImageClose (outImage, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, outImage->name);
    return;
  }
  /* Free image buffer if not memory resident  */
  if (outImage->mySel->FileType!=OBIT_IO_MEM) 
    outImage->image = ObitFArrayUnref(outImage->image);

} /* end  ObitImageUtilPBApply */

/**
 * Make an image of the antenna primary beam pattern based on the pointing
 * position in an image.
 * For frequencies < 1 GHz uses the VLA polynomial gain curves,
 * for higher frequencies, it uses a jinc function based on the antenna size.
 * \param pntImage Image with pointing position
 * \param outImage Image to be written.  Must be previously instantiated.
 * \param outPlane  Desired plane in outImage on planes 3-5; ignored if memOnly
 * \param antSize  Antenna size
 * \param minGain  Min. allowed antenna gain, lower values are blanked
 * \param err      Error stack, returns if not empty.
 */
void 
ObitImageUtilPBImage (ObitImage *pntImage, ObitImage *outImage, 
		     olong *outPlane, ofloat antSize, ofloat minGain, ObitErr *err)
{
  ObitIOSize IOBy;
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong ix, iy, indx, pos[2];
  ofloat inPixel[2], *out, fblank = ObitMagicF();
  odouble RAPnt, DecPnt, Freq, ra, dec, xx, yy, zz, dist ;
  ofloat pbf, equinox;
  gboolean doJinc, bad;
  ObitImageDesc *outDesc;
  gchar *routine = "ObitImageUtilPBApply";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(pntImage));
  g_assert (ObitImageIsA(outImage));

  if (antSize<0.01) antSize = 25.0; /* default antenna size */

  /* Get pointing position */
  pntImage->extBuffer = TRUE;  /* Don't need buffer */
  if ((ObitImageOpen (pntImage, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, pntImage->name);
    return;
  }

  /* Use  "Observed" position if given */
  RAPnt   = pntImage->myDesc->obsra;
  DecPnt  = pntImage->myDesc->obsdec;
  equinox = pntImage->myDesc->equinox ;
  if ((abs(RAPnt)<1.0e-5) && (abs(DecPnt)<1.0e-5)) {
    /* if zeroes - use reference position */
    RAPnt  = pntImage->myDesc->crval[pntImage->myDesc->jlocr];
    DecPnt = pntImage->myDesc->crval[pntImage->myDesc->jlocd];
  }
  RAPnt  *= DG2RAD;
  DecPnt *= DG2RAD;
  if ((ObitImageClose (pntImage, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, pntImage->name);
    return;
  }
  pntImage->extBuffer = FALSE;  /* May need buffer later */

  /* Do I/O by plane and all of plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListPut (outImage->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
  dim[0] = 7;
  for (i=0; i<5; i++) blc[i+2] = trc[i+2] = outPlane[i];
  ObitInfoListPut (outImage->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (outImage->info, "TRC", OBIT_long, dim, trc, err);

  /* Open image */
  if ((ObitImageOpen (outImage, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, outImage->name);
    return;
  }
  
  /* check coordinate types */
  if (outImage->myDesc->coordType != pntImage->myDesc->coordType) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Output (%s) and pointing (%s) images have different coordinate types ", 
		   routine, outImage->name, pntImage->name);
    return;
  }

  /* Precess pointing position is needed */
  if ((abs(equinox-1950.0)<0.01) && 
      (abs(outImage->myDesc->equinox-2000.0)<0.01))
    ObitSkyGeomBtoJ (&RAPnt, &DecPnt);
  else if ((abs(equinox-2000.0)<0.01) && 
      (abs(outImage->myDesc->equinox-1950.0)<0.01))
    ObitSkyGeomJtoB (&RAPnt, &DecPnt);

  /* Get output aray pointer */
  pos[0] = pos[1] = 0;
  out = ObitFArrayIndex (outImage->image, pos);

  /* Set up - get frequency, default 1 GHz*/
  Freq = 1.0e9;
  if (outImage->myDesc->jlocf>=0) 
    Freq = outImage->myDesc->crval[outImage->myDesc->jlocf];
  /* which beam model to use */
  doJinc = (Freq >= 1.0e9);
  outDesc = outImage->myDesc; /* output descriptor */

  /* Loop over image  */
  for (iy = 1; iy<=outDesc->inaxes[1]; iy++) { /* loop in y */
    inPixel[1] = (ofloat)iy;
    for (ix = 1; ix<=outDesc->inaxes[0]; ix++) {/* loop in x */
      inPixel[0] = (ofloat)ix;

      /* array index in in and out for this pixel */
      indx = (iy-1) * outDesc->inaxes[0] + (ix-1);

      /* Convert pixel to position */
      bad = 
	ObitSkyGeomWorldPos(inPixel[0], inPixel[1],
			    outDesc->crval[outDesc->jlocr], outDesc->crval[outDesc->jlocd],
			    outDesc->crpix[outDesc->jlocr], outDesc->crpix[outDesc->jlocd],
			    outDesc->cdelt[outDesc->jlocr], outDesc->cdelt[outDesc->jlocd],
			    outDesc->crota[outDesc->jlocd], &outDesc->ctype[outDesc->jlocr][4],
			    &ra, &dec);
      if (bad!=0) {
	Obit_log_error(err, OBIT_Error, 
		       "%s: Error %d determining location of pixel in %s", 
		       routine, bad, outImage->name);
	return;
      }
      
      /* Separation from pointing center */
      xx = DG2RAD * ra;
      yy = DG2RAD * dec;
      zz = sin (yy) * sin (DecPnt) + cos (yy) * cos (DecPnt) * cos (xx-RAPnt);
      zz = MIN (zz, 1.000);
      dist = acos (zz) * RAD2DG;
      
      /* primary beam correction */
      if (doJinc) {
	pbf = ObitPBUtilJinc (dist, Freq, antSize, 0.0);
      } else {
	pbf = ObitPBUtilPoly (dist, Freq, 0.0);
      } 
      
      /* debug
	 if (pbf>0.9) {
	 fprintf (stderr,"pdf %f @ %d  %d %lf\n",pbf,ix,iy,dist* RAD2DG*3600.0);
	 } */

      /* Clip by minGain */
      if (pbf<minGain) pbf = fblank;
      
      /* Save beam image */
      out[indx] = pbf;

    } /* end loop over x */
  } /* end loop over y */
  

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
  /* Free image buffer  if not memory resident */
  if (outImage->mySel->FileType!=OBIT_IO_MEM) 
    outImage->image = ObitFArrayUnref(outImage->image);

} /* end  ObitImageUtilPBImage */

/**
 * Use maximum baseline length and maximum W to set imaging cell size
 * and maximum undistorted field of view
 *  If Cells is specified, then it is used, otherwise it is determined from the 
 *  longest baseline in the uvdata (MaxBL)
 *  If input value of Radius is given, it is used as the field size.  
 *  Otherwise, maximum field size from Lecture 2  (A. R. Thompson) in 
 *  "Synthesis Imaging in Radio Astronomy II",  PASP, vol. 180, 1999, p. 24 is used.
 *  A correction is applied for the actual range in W.
 * \param MaxBL   maximum baseline length (sqrt(u*u+v*v))
 * \param MaxW    Max abs(w) in data.
 * \param Cells   Cell spacing in asec.  If zero on input the value is set
 *                based on MaxBL (1/4 min fringe spacing)
 * \param Radius  Maximum undistorted Field of view in cells (Cells).
 *                If zero on input, the value is set from MaxW.
 */
void ObitImageUtilImagParm (ofloat MaxBL, ofloat MaxW,
			    ofloat *Cells, ofloat *Radius)
{
  ofloat fs, hpbw, maxfr;
  
  /* Set values if not defined */
  /* Cell spacing based on MAXBL */
  if (*Cells<=0.0) {
    /* Maximum fringe spacing in asec - fudge a bit to avoid weighting problem */
    fs = RAD2AS / (MaxBL*1.0001);
    *Cells = fs / 4.0;
  }
  
  /* Field size limited by W  distortion  */
  if (*Radius<=0.0) {
    /* Estimate beam size */
    hpbw = (4.0 * (*Cells)) * AS2RAD;
    /* Maximum field size in radians */
    maxfr = 0.33 * sqrt (hpbw);
    /* Correct by sqrt ratio of MaxBL to  MaxW */
    /*DEBUGmaxfr = maxfr * MaxBL / MaxW; */
    maxfr = maxfr * sqrt(MaxBL / MaxW);
    /*  Assume circular clean beam - Undisturbed field size */
    *Radius =  maxfr*RAD2AS / (*Cells) + 0.5;
    /*DEBUG *Radius =  (0.5 * maxfr*RAD2AS) / (*Cells) + 0.5;*/
  }

} /* end ObitImageUtilImagParm */

/**
 * Rudimentry ObitFArray to ObitImage converter
 * Create FITS image file with name and disk given,
 * Makes basic image header from inArray and writes the file
 * as a floating image.
 * Currently only does two to four dimensions.
 * \param fileName output FITS image file name
 * \param disk     output FITS disk number
 * \param inArray  Data array to be written
 * \param err      Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitImage* ObitImageUtilArray2Image (gchar *fileName, olong disk, 
				     ObitFArray *inArray, ObitErr *err)
{
  ObitImage *out=NULL;
  ObitImageDesc *myDesc;
  ofloat *data;
  olong pos[5]={0,0,0,0,0};
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong i, j, planeNo[5] = {1,1,1,1,1};
  gchar *routine = "ObitImageUtilArray2Image";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitFArrayIsA(inArray));

  /* Create basic output image */
  out = newObitImage(fileName);

  /* Define file */
  ObitImageSetFITS(out,OBIT_IO_byPlane,disk,fileName,blc,trc,err);
  if (err->error) Obit_traceback_val (err, routine, fileName, out);

  /* Create rudimentary header */
  myDesc = out->myDesc;
  myDesc->bitpix = -32;  /* Floating */
  myDesc->naxis = inArray->ndim;
  g_snprintf (myDesc->object, IMLEN_VALUE-1, "%s", fileName);
  for (i=0; i<myDesc->naxis; i++) { /* Axis descriptors */
    myDesc->inaxes[i] = inArray->naxis[i];
    myDesc->crpix[i] = 1.0 + inArray->naxis[i]/2.0;
    myDesc->crota[i] = 0.0;
    myDesc->cdelt[i] = 1.0/3600.0;
    myDesc->crval[i] = 0.0;
    g_snprintf (myDesc->ctype[i], IMLEN_KEYWORD-1, "AXIS%d",i+1);
  }

  /* Open file */
  ObitImageOpen (out, OBIT_IO_WriteOnly, err);
  if (err->error) Obit_traceback_val (err, routine, fileName, out);

  /* Loop over fourth dimension writing planes */
  for (j=0; j<MAX(1,inArray->naxis[3]); j++) {
    /* Loop over third */
    for (i=0; i<MAX(1,inArray->naxis[2]); i++) {
      /* pointer to start of plane*/
      pos[0] = 0; pos[1] = 0; pos[2] = i; pos[3] = j;
      data = ObitFArrayIndex(inArray, pos); 
      /* Write it */
      planeNo[0] = i+1; planeNo[1] = j+1; 
      ObitImagePutPlane (out, data, planeNo, err);
      if (err->error) Obit_traceback_val (err, routine, fileName, out);
    } /* end loop over third axis */
  } /* end loop writing planes */

  /* Close file */
  ObitImageClose (out, err);
  if (err->error) Obit_traceback_val (err, routine, fileName, out);
  /* Free image buffer if not memory resident */
  if (out->mySel->FileType!=OBIT_IO_MEM) 
    out->image = ObitFArrayUnref(out->image);

  /* Give message */
  Obit_log_error(err, OBIT_InfoErr, 
		 "%s: Wrote ObitFArray %s to FITS image %s disk %d", 
		 routine, inArray->name, fileName, disk);
  return out;
} /* end  ObitImageUtilArray2Image */

/**
 * Quantize an image  at a specified quantization level or fraction of 
 * the image RMS and write to a FITS image.
 * Image RMS derived from histogram fitting and should be a reasonable 
 * estimate of the "noise".
 * Selection by blc, trc in inImage is honored.
 * \param inImage  Image to quantize, parameters in info:
 * \li "factor" OBIT_float (1,1,1) quantize at factor*RMS [def 0.2]
 *              RMS is the minimum rms in any selected plane.
 * \li "quant"  OBIT_float (1,1,1) quantization level, 
 *              has presidence over factor, def.(or <=0) use factor
 * \param fileName output FITS image file name
 * \param disk     output FITS directory number
 * \param err      Error stack, returns if not empty.
 * \return pointer to the new object, may be NULL on failure.
 */
ObitImage* ObitImageUtilQuanFITS (ObitImage *inImage, gchar *fileName, 
				  olong disk, ObitErr *err)
{
  ObitImage *out=NULL;
  ObitIOSize IOBy;
  ObitImageDesc *outDesc;
  ObitIOCode iretCode, oretCode;
  ofloat factor, quant;
  ofloat planeRMS, minRMS, planeMin, minMin, planeMax, maxMax, dr; 
  ObitInfoType type;
  gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong pos[5]={0,0,0,0,0};
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong bitpix, tt, bb;
  gchar *today=NULL;
  gchar *routine = "ObitImageUtilQuanFITS";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitImageIsA(inImage));

  /* Control parameters */
  factor = 0.2;
  ObitInfoListGetTest(inImage->info, "factor", &type, dim, &factor);
  quant = -1.0;
  ObitInfoListGetTest(inImage->info, "quant", &type, dim, &quant);

  /* DEBUG
     fprintf (stderr, "DEBUG %s factor %f quant %f\n",routine, factor, quant); */

  /* Do I/O by plane (but keep any blc, trc) */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListPut (inImage->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
  inImage->extBuffer = FALSE;  /* Make sure it has buffer */

  /* Need statistics? */
  if (quant<=0.0) {
    /* DEBUG 
    fprintf (stderr, "DEBUG %s getting statistics\n",routine);*/

    /* Open input image */
    iretCode = ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
    if (err->error) Obit_traceback_val (err, routine, inImage->name, out);
    
    /* Loop through input image finding minimum plane RMS, max, min */
    minRMS = 1.0e25;
    minMin = 1.0e25;
    maxMax =-1.0e25;
    while (iretCode==OBIT_IO_OK) {
      iretCode = ObitImageRead (inImage, NULL, err);
      if (iretCode == OBIT_IO_EOF) break;  /* Done */
      if (err->error) Obit_traceback_val (err, routine, inImage->name, out);
      
      /* Get plane statistics */
      planeRMS = ObitFArrayRMS (inImage->image);
      planeMin = ObitFArrayMin (inImage->image, pos);
      planeMax = ObitFArrayMax (inImage->image, pos);
      if (fabs(planeRMS)>0.0) {  /* Ignore empty planes */
	minRMS = MIN (minRMS, planeRMS);
	minMin = MIN (minMin, planeMin);
	maxMax = MAX (maxMax, planeMax);
      }
    } /* end loop reading */
    
    iretCode = ObitImageClose (inImage, err);  /* Close input */
    if (err->error) Obit_traceback_val (err, routine, inImage->name, out);

    /* Set quantization level */
    /* Deal with zero RMS (probably poor statistics) - use 1/32700 peak */
    if (minRMS>0.0) 
      quant = factor * minRMS;
    else
      quant =  MAX (fabs(maxMax), fabs(minMin)) / 32700.0;
    /* end get statistics and quantization level */
  } else {
    /* Get max, min from header */
    minMin = inImage->myDesc->minval;
    maxMax = inImage->myDesc->maxval;
  }

  /* Set bitpix based on the dynamic range needed */
  dr = MAX (fabs(maxMax), fabs(minMin)) / quant;
  /* DEBUG 
     fprintf (stderr, "DEBUG %s dynamic range %f range %f %f\n",routine, dr,maxMax, minMin );*/
  if (dr < 32760.0)  bitpix = 16;
  else if (dr < 2147483600.0) bitpix = 32;
  else { /* Can't write as integer */
    Obit_log_error(err, OBIT_Error, 
		   "%s: cannot write %s with quantization %f as integer", 
		   routine, fileName, quant);
    return out;
  }

  /* Tell what's going on */
  Obit_log_error(err, OBIT_InfoErr, 
		 "Writing %s to FITS image %s disk %d", 
		 inImage->name, fileName, disk);
  Obit_log_error(err, OBIT_InfoErr, 
		 "quantization=%f bitpix=%d", quant, bitpix);

  /* Create basic output image */
  out = newObitImage(fileName);

  /* Define file */
  ObitImageSetFITS(out,OBIT_IO_byPlane,disk,fileName,blc,trc,err);
  if (err->error) Obit_traceback_val (err, routine, fileName, out);

  /* Copy header descriptive material */
  ObitImageDescCopyDesc (inImage->myDesc, out->myDesc, err);
  if (err->error) Obit_traceback_val (err, routine, inImage->name, out);
  outDesc = out->myDesc;
  outDesc->bitpix = bitpix;
  outDesc->minval = minMin;
  outDesc->maxval = maxMax;
  g_snprintf (outDesc->origin, IMLEN_VALUE-1, "Generated by Obit");
  /* Save quantization value */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(outDesc->info, "Quant", OBIT_float, dim, &quant);

  /* Creation date today */
  today = ObitToday();
  strncpy (outDesc->date, today, IMLEN_VALUE-1);
  if (today) g_free(today);
 
  /* Set size of output */
  outDesc->naxis = inImage->myDesc->naxis;
  /* Correct for selection */
  for (i=0; i<outDesc->naxis; i++) { 
    tt = inImage->mySel->trc[i];
    if (tt<=0) tt = inImage->myDesc->inaxes[i];
    bb = inImage->mySel->blc[i];
    if (bb<=0) bb = 1;
    outDesc->inaxes[i] = tt-bb+1;
    outDesc->crpix[i]  = inImage->myDesc->crpix[i] - bb + 1;
  }

  /* Fully define */
  out->extBuffer = TRUE;  /* Don't need to assign buffer here */
  /* Open and close */
  ObitImageOpen(out, OBIT_IO_WriteOnly, err);
  /* Force scaling */
  out->myStatus = OBIT_Modified;

  ((ObitIOImageFITS*)out->myIO)->dataMod = TRUE;
  /* update scaling */
  ObitIOImageFITSUpdateScale (((ObitIOImageFITS*)out->myIO),quant,  err);
  /*   ObitImageClose(out, err);*/
  if (err->error) Obit_traceback_val (err, routine, fileName, out);
  
  /* Use external buffer for writing output */
  out->extBuffer = TRUE;

  /* Open files */
  iretCode = ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_val (err, routine, inImage->name, out);
  /* oretCode = ObitImageOpen (out, OBIT_IO_WriteOnly, err); */
  if (err->error) Obit_traceback_val (err, routine, fileName, out);

  /* Loop copying planes */
  oretCode = OBIT_IO_OK;
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    iretCode = ObitImageRead (inImage, NULL, err);
    if (iretCode == OBIT_IO_EOF) break;  /* Done */
    if (err->error) Obit_traceback_val (err, routine, inImage->name, out);
    /* Write plane */
    oretCode = ObitImageWrite(out, inImage->image->array, err);
    if (err->error) Obit_traceback_val (err, routine, fileName, out);
  } /* end loop writing planes */

  /* Close files */
  iretCode = ObitImageClose (inImage, err);  /* Close input */
  if (err->error) Obit_traceback_val (err, routine, inImage->name, out);
  oretCode = ObitImageClose (out, err);
  if (err->error) Obit_traceback_val (err, routine, fileName, out);
  /* Unset external buffer for writting */
  out->extBuffer = FALSE;

  /* Release Input image buffer  if not memory resident */
  if (inImage->mySel->FileType!=OBIT_IO_MEM) 
    inImage->image = ObitFArrayUnref(inImage->image);

  return out;
} /* end  ObitImageUtilQuanFITS */


/**
 * Create an image and fill the descriptor values for an image cube
 * based on  the descriptor for a single plane and for the uv data 
 * creating the image.
 * This should be called before the image is Opened or instantiated.
 * \param inDesc    Input Image Descriptor.
 * \param UVDesc    Input UV Descriptor.
 * \param outDesc   Output Image Descriptor 
 * \param Stokes    Stokes parameter of image ' '=>'I', (I, Q, U, V, R, L)
 * \param bchan     first (1-rel) channel in UVDesc
 * \param echan     highest (1-rel) channel in UVDesc
 * \param incr      channel increment in input
 * \param err       Error stack, returns if not empty.
 */
void 
ObitImageUtilMakeCube (ObitImageDesc *inDesc, ObitUVDesc *UVDesc, 
		       ObitImageDesc *outDesc, 
		       gchar *Stokes, olong bchan, olong echan, olong incr, ObitErr *err)
{
  olong numberChann;
  gchar *name, *today=NULL;
  gchar *routine = "ObitImageUtilMakeCube";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageDescIsA(inDesc));
  g_assert (ObitUVDescIsA(UVDesc));

  /* Save output name */
  if (outDesc->name) name = g_strdup (outDesc->name);
  else  name = g_strdup ("Descriptor");

  /* Most info from inDesc */
  outDesc = ObitImageDescCopy (inDesc, outDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, inDesc->name);

  /* Creation date today */
  today = ObitToday();
  strncpy (outDesc->date, today, IMLEN_VALUE-1);
  if (today) g_free(today);
 
  /* restore name */
  if (outDesc->name) g_free(outDesc->name);
  outDesc->name = name;

  /* Set number of channels - include effects of averaging*/
  numberChann = MIN (echan, UVDesc->inaxes[UVDesc->jlocf]) - MAX (1, bchan) + 1;
  numberChann = (olong)((((ofloat)numberChann) / MAX (1, incr)) + 0.999);
  outDesc->inaxes[outDesc->jlocf] = MAX (1, numberChann);

  /* Stokes parameter */
  if ((Stokes[0]=='I') || (Stokes[0]==' ')) outDesc->crval[outDesc->jlocs] = 1.0;
  else if (Stokes[0]=='Q') outDesc->crval[outDesc->jlocs] =  2.0;
  else if (Stokes[0]=='U') outDesc->crval[outDesc->jlocs] =  3.0;
  else if (Stokes[0]=='V') outDesc->crval[outDesc->jlocs] =  4.0;
  else if (Stokes[0]=='R') outDesc->crval[outDesc->jlocs] = -1.0;
  else if (Stokes[0]=='L') outDesc->crval[outDesc->jlocs] = -1.0;

  /* reset image max/min */
  outDesc->maxval    = -1.0e20;
  outDesc->minval    =  1.0e20;

  return;
} /* end ObitImageUtilMakeCube */

/**
 * Write the (first) plane from image in to a plane in out.
 * \param in        Input image with plane to copy
 * \param out       Output cube to accept plane
 * \param plane     (1-rel) pixel indices for planes 3-7 in out.
 * \param err       Error stack, returns if not empty.
 */
void ObitImageUtilInsertPlane (ObitImage *in, ObitImage *out, olong *plane, 
			       ObitErr *err)
{
  ObitIOSize IOBy = OBIT_IO_byPlane;
  ObitImageDesc *inDesc=NULL, *outDesc=NULL;
  olong i;
  olong blc[IM_MAXDIM] = {1,1,1,1,1};
  olong trc[IM_MAXDIM] = {1,1,1,1,1};
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "ObitImageUtilInsertPlane";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageIsA(in));
  g_assert (ObitImageIsA(out));
  g_assert (plane!=NULL);

  inDesc = in->myDesc;
  outDesc = out->myDesc;
  /* Check size of planes */
  Obit_return_if_fail(((inDesc->inaxes[0]==outDesc->inaxes[0]) && 
		       (inDesc->inaxes[1]==outDesc->inaxes[1])), err,
		      "%s: Image planes incompatible  %d!= %d or  %d!= %d", 
		      routine, inDesc->inaxes[0], outDesc->inaxes[0], 
		      inDesc->inaxes[1], outDesc->inaxes[1]) ;
 
  Obit_return_if_fail(((plane[0]<=outDesc->inaxes[2]) && 
		       (plane[1]<=outDesc->inaxes[3])), err,
		      "%s: Output does not have plane %d  %d", 
		      routine, plane[0], plane[1]);

  /* Read input plane */
  /* Set blc, trc */
  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<2; i++) trc[i] = inDesc->inaxes[i];
  for (i=2; i<IM_MAXDIM; i++) trc[i] = 1;
  
  dim[0] = 1;
  ObitInfoListPut (in->info, "IOBy", OBIT_long, dim, &IOBy, err);
  dim[0] = 7;
  ObitInfoListPut (in->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (in->info, "TRC", OBIT_long, dim, trc, err);
  ObitImageOpen (in, OBIT_IO_ReadOnly, err);
  ObitImageRead (in, NULL, err);
  ObitImageClose (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Write to output */
  /* Set blc, trc */
  for (i=0; i<2; i++) blc[i] = 1;
  for (i=2; i<IM_MAXDIM; i++) blc[i] = MAX (1,plane[i-2]);
  for (i=0; i<2; i++) trc[i] = outDesc->inaxes[i];
  for (i=2; i<IM_MAXDIM; i++) trc[i] = MAX (1,plane[i-2]);
  dim[0] = 1;
  ObitInfoListPut (out->info, "IOBy", OBIT_long, dim, &IOBy, err);
  dim[0] = 7;
  ObitInfoListPut (out->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (out->info, "TRC", OBIT_long, dim, trc, err);
  out->extBuffer = TRUE;  /* Don't need output buffer */
  ObitImageOpen (out, OBIT_IO_ReadWrite, err);
  ObitImageWrite (out, in->image->array, err);
  /* Copy Beam information if needed */
  if ((out->myDesc->beamMaj<=0.0) || (out->myDesc->beamMin<=0.0)) {
    out->myDesc->beamMaj = in->myDesc->beamMaj;
    out->myDesc->beamMin = in->myDesc->beamMin;
    out->myDesc->beamPA  = in->myDesc->beamPA;
    out->myDesc->niter = 1;
  }
  ObitImageClose (out, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);
  out->extBuffer = FALSE;  /* May need buffer later */

  /* free image memory if not memory resident */
  if (in->mySel->FileType!=OBIT_IO_MEM) 
    in->image = ObitFArrayUnref(in->image);

} /* end ObitImageUtilInsertPlane */

/**
 * Insert multiple planes from image in starting at plane in out.
 * \param in        Input image cube with planes to copy
 *                  Any BLC, TRC are honored
 * \param out       Output cube to accept planes
 * \param plane     (1-rel) pixel indices for planes 3-7 in out.
 * \param axExp     (1-rel) axis number being expanded (usually 3)
 * \param err       Error stack, returns if not empty.
 */
void ObitImageUtilInsertCube (ObitImage *in, ObitImage *out, olong *plane, 
			      olong axExp, ObitErr *err)
{
  ObitIOSize IOBy = OBIT_IO_byPlane;
  ObitImageDesc *inDesc=NULL, *outDesc=NULL;
  olong i, iplane, nplane;
  olong inblc[IM_MAXDIM]={1,1,1,1,1,1,1}, blc[IM_MAXDIM]={1,1,1,1,1,1,1};
  olong intrc[IM_MAXDIM]={0,0,0,0,0,0,0}, trc[IM_MAXDIM]={0,0,0,0,0,0,0};
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gchar *routine = "ObitImageUtilInsertPlane";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageIsA(in));
  g_assert (ObitImageIsA(out));
  g_assert (plane!=NULL);

  /* Access images a plane at a time */
  dim[0] = 1;
  ObitInfoListPut (in->info, "IOBy", OBIT_long, dim, &IOBy, err);
  ObitInfoListPut (out->info, "IOBy", OBIT_long, dim, &IOBy, err);

  inDesc  = in->myDesc;
  outDesc = out->myDesc;

  /* Open input to ensure size set OK */
  ObitImageOpen (in, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Get selected input size */
  ObitInfoListGetTest (in->info, "BLC", &type, dim, inblc); 
  ObitInfoListGetTest (in->info, "TRC", &type, dim, intrc);
  for (i=0; i<IM_MAXDIM; i++) {
    inblc[i] = MAX(1, MIN(inblc[i],inDesc->inaxes[i]));
    if (intrc[i]<=0) intrc[i] = inDesc->inaxes[i];
    intrc[i] = MAX(1, MIN(intrc[i],inDesc->inaxes[i]));
  }

  /* Check size of planes */
  Obit_return_if_fail(((inDesc->inaxes[0]==outDesc->inaxes[0]) && 
		       (inDesc->inaxes[1]==outDesc->inaxes[1])), err,
		      "%s: Image planes incompatible  %d!= %d or  %d!= %d", 
		      routine, inDesc->inaxes[0], outDesc->inaxes[0], 
		      inDesc->inaxes[1], outDesc->inaxes[1]) ;
 
  Obit_return_if_fail(((axExp>2) && (axExp<=IM_MAXDIM)), err,
		      "%s: Illegal axis to expand %d",  routine, axExp);
 
  Obit_return_if_fail(((plane[0]<=outDesc->inaxes[2]) && 
		       (plane[1]<=outDesc->inaxes[3])), err,
		      "%s: Output does not have plane %d  %d", 
		      routine, plane[0], plane[1]);

  /* How many input planes? */
  nplane = inDesc->inaxes[2];

  /* Loop over output planes */
  for (iplane=0; iplane<nplane; iplane++) {

    /* Read input plane */
    ObitImageRead (in, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Open output - Set blc, trc */
    for (i=0; i<2; i++) blc[i] = 1;
    for (i=2; i<IM_MAXDIM; i++) blc[i] = MAX (1,plane[i-2]);
    blc[axExp-1] += iplane;
    for (i=0; i<2; i++) trc[i] = outDesc->inaxes[i];
    for (i=2; i<IM_MAXDIM; i++) trc[i] = MAX (1,plane[i-2]);
    trc[axExp-1] += iplane;
    dim[0] = 7;
    ObitInfoListPut (out->info, "BLC", OBIT_long, dim, blc, err); 
    ObitInfoListPut (out->info, "TRC", OBIT_long, dim, trc, err);
    out->extBuffer = TRUE;  /* Don't need output buffer */
    ObitImageOpen (out, OBIT_IO_ReadWrite, err);
  
   /* Write to output */
    ObitImageWrite (out, in->image->array, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    /* Copy Beam information if needed */
    if ((out->myDesc->beamMaj<=0.0) || (out->myDesc->beamMin<=0.0)) {
      out->myDesc->beamMaj = in->myDesc->beamMaj;
      out->myDesc->beamMin = in->myDesc->beamMin;
      out->myDesc->beamPA  = in->myDesc->beamPA;
      out->myDesc->niter = 1;
    }
    ObitImageClose (out, err);
    if (err->error) Obit_traceback_msg (err, routine, out->name);
    out->extBuffer = FALSE;  /* May need buffer later */
  } /* end loop over planes */

  /* Close input */
  ObitImageClose (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Reset input BLC, TRC */
  dim[0] = 7;
  ObitInfoListAlwaysPut (in->info, "BLC", OBIT_long, dim, inblc); 
  ObitInfoListAlwaysPut (in->info, "TRC", OBIT_long, dim, intrc);

  /* free image memory if not memory resident */
  if (in->mySel->FileType!=OBIT_IO_MEM) 
    in->image = ObitFArrayUnref(in->image);

} /* end ObitImageUtilInsertCube */

/**
 * Fill in an image Descriptor from a UV Descriptor.
 * Needs any xshift an y shift filled into the image prior to call.
 * Information about the first two axes other than the type an 
 * coordinate value need to be set separately.
 * to get the final position correct.
 * \param UVDesc    Input UV Descriptor.
 * \param imageDesc Output image Descriptor
 * \param nchavg    How many uv channels to average per image channel.
 *                  Ignored if uv data has multiple IFs.
 */
void 
ObitImageUtilUV2ImageDesc(ObitUVDesc *UVDesc, ObitImageDesc *imageDesc, 
			  olong nchavg)
{
  olong i, iaxis, nif, nfreq, nch;
  odouble sum;
  gchar *st1, *st2;
  gchar *today=NULL;

  /* error checks */
  g_assert (ObitUVDescIsA(UVDesc));
  g_assert (ObitImageDescIsA(imageDesc));
  
  /* Be sure UV descriptor is indexed */
  ObitUVDescIndex(UVDesc);

  /* Creation date today */
  today = ObitToday();
  strncpy (imageDesc->date, today, IMLEN_VALUE-1);
  if (today) g_free(today);
 
  /* loop over axes */

  iaxis = 0;
  /* RA axis, inaxes, cdelt, crota, xshift set else where */
  /* Form label string */
  st1 = imageDesc->ctype[iaxis];
  st2 = UVDesc->ctype[UVDesc->jlocr];
  for (i=0; i<4; i++) { /* Axis type */
    st1[i]=st2[i];
    if (st1[i]==' ') st1[i]='-';
    if (st1[i]==0)   st1[i]='-';
  }
  st2 = UVDesc->ptype[UVDesc->ilocu]; /* Projection */
  for (i=0; i<4; i++)  st1[i+4]=st2[i+4];
  st1[9] = 0;

  /* Reference pixel */
  imageDesc->crpix[iaxis] = 1.0 + imageDesc->inaxes[iaxis] / 2.0;

  /* Dec axis, inaxes, cdelt, crota, xshift set else where */
  iaxis++;
  /* Form label string */
  st1 = imageDesc->ctype[iaxis];
  st2 = UVDesc->ctype[UVDesc->jlocd];
  for (i=0; i<4; i++) { /* Axis type */
    st1[i]=st2[i];
    if (st1[i]==' ') st1[i]='-';
    if (st1[i]==0)   st1[i]='-';
  }
  st2 = UVDesc->ptype[UVDesc->ilocu]; /* Projection */
  for (i=0; i<4; i++)  st1[i+4]=st2[i+4];
  st1[9] = 0;

  /* Reference pixel */
  imageDesc->crpix[iaxis] = 1.0 + imageDesc->inaxes[iaxis] / 2.0;

  /* Position from uv reference and shift */
  ObitSkyGeomXYShift (UVDesc->crval[UVDesc->jlocr], 
		      UVDesc->crval[UVDesc->jlocd],
		      UVDesc->xshift, UVDesc->yshift, 
		      imageDesc->crota[1],
		      &imageDesc->crval[0], &imageDesc->crval[1]);

  /* Frequency Axis */
  iaxis++;
  /* How many? */
  if (UVDesc->jlocif>=0) nif = UVDesc->inaxes[UVDesc->jlocif];
  else  nif = 1;
  nfreq = UVDesc->inaxes[UVDesc->jlocf];
  /* Initially set for continuum */
  strncpy (imageDesc->ctype[iaxis], "FREQ    ", IMLEN_KEYWORD-1);
  imageDesc->inaxes[iaxis] = 1;  /* Only one for continuum */
  imageDesc->crpix[iaxis] = 1.0; /* reference pixel */
  imageDesc->crota[iaxis] = 0.0; /* no possible meaning */
  /* coordinate increment = total bandwidth */
  imageDesc->cdelt[iaxis] = nif*nfreq*UVDesc->cdelt[UVDesc->jlocf]; 
  /* Output frequency is average frequency */
  sum = 0.0;
  for (i=0; i<nif*nfreq; i++) sum += UVDesc->freqArr[i];
  imageDesc->crval[iaxis] = sum / (nif*nfreq);

  /* More complex if spectral line observations */
  /* May average channels */
  if ((nif==1) && (nchavg<nfreq)) {
    nch = nfreq / nchavg; /* how many output channels? */
    imageDesc->inaxes[iaxis] = nch;
    imageDesc->cdelt[iaxis] = nch*UVDesc->cdelt[UVDesc->jlocf]; 
    imageDesc->crpix[iaxis] = 1.0 + (nch-1.0) / 2.0; /* reference pixel */
    /* average first nch uv channels out output reference frequency */
    sum = 0.0;
    for (i=0; i<nch; i++) sum += UVDesc->freqArr[i];
    imageDesc->crval[iaxis] = sum / (nch);
} /* end update for spectral line */

  /* Stokes Axis */
  iaxis++;
  imageDesc->inaxes[iaxis] = 1;  /* Only one */
  strncpy (imageDesc->ctype[iaxis], "STOKES  ", IMLEN_KEYWORD-1);
  imageDesc->crval[iaxis] = UVDesc->crval[UVDesc->jlocs];
  imageDesc->crpix[iaxis] = 1.0; /* reference pixel */
  imageDesc->cdelt[iaxis] = 1.0; /* coordinate increment */
  imageDesc->crota[iaxis] = 0.0; /* no possible meaning */

  /* Total number of axes */
  imageDesc->naxis = iaxis+1;

  /* Copy information not directly related to an axis */
  /* Strings */
  strncpy (imageDesc->object, UVDesc->object, IMLEN_VALUE-1);
  strncpy (imageDesc->teles,  UVDesc->teles,  IMLEN_VALUE-1);
  strncpy (imageDesc->instrument,  UVDesc->instrument,  IMLEN_VALUE-1);
  strncpy (imageDesc->observer,  UVDesc->observer,  IMLEN_VALUE-1);
  strncpy (imageDesc->obsdat, UVDesc->obsdat, IMLEN_VALUE-1);
  strncpy (imageDesc->origin, UVDesc->origin, IMLEN_VALUE-1);
  strncpy (imageDesc->bunit,  "JY/BEAM ",     IMLEN_VALUE-1);
  /* Set current date */
  ObitImageUtilCurDate (imageDesc->date, IMLEN_VALUE-1);

  /* Other */
  imageDesc->altCrpix     = UVDesc->altCrpix;
  imageDesc->altRef       = UVDesc->altRef;
  imageDesc->restFreq     = UVDesc->restFreq;
  imageDesc->VelReference = UVDesc->VelReference;
  imageDesc->VelDef       = UVDesc->VelDef;
  imageDesc->xshift       = UVDesc->xshift;
  imageDesc->yshift       = UVDesc->yshift;
  imageDesc->epoch        = UVDesc->epoch;
  imageDesc->equinox      = UVDesc->equinox;
  imageDesc->obsra        = UVDesc->obsra;
  imageDesc->obsdec       = UVDesc->obsdec;

  /* initialize some values */
  imageDesc->areBlanks = FALSE;
  imageDesc->niter     = 0;
  imageDesc->maxval    = -1.0e20;
  imageDesc->minval    =  1.0e20;
  imageDesc->bitpix    = -32;
  imageDesc->beamMaj   = 0.0;
  imageDesc->beamMin   = 0.0;
  imageDesc->beamPA    = 0.0;

  /* Index Image descriptor */
  ObitImageDescIndex(imageDesc);

} /* end ObitImageUtilUV2ImageDesc */

/**
 * Make flux weighted velocity image from a velocity cube.
 * Input image is clipped to only significant pixels.
 * Convolution of each plane by Parms[3] cells used to mask image
 * \param inImage  Input velocity cube image
 * Parameters in info:
 * \li "BLC"     OBIT_int (7) Lowest x,y,v pixel number selected [def 1,1,1]
 * \li "TRC"     OBIT_int (7) Highest x,y,v pixel number selected [def all]
 *      Note: the first two output axes will have an even number of pixels.
 * \li "Parms"   OBIT_float (4) Parameters
 *     [0] min. RMS (convolved image)
 *     [1] min. fraction of  peak
 *     [2] <0.5 => blank fill, else zero fill
 *     [3] Convolution size (0->5)
 * \param outImage Image to be written.  Must be previously instantiated.
 * \param err      Error stack, returns if not empty.
 */
void 
ObitImageUtilVel (ObitImage *inImage, ObitImage *outImage, ObitErr *err)
{
  ObitIOCode   iretCode, oretCode;
  ofloat       RMS, maxF, newVal, *Parms=NULL, fblank =  ObitMagicF();
  ofloat       vel, minAllow, FWHM;
  olong         iplane;
  ObitFArray   *accum=NULL, *accumWt=NULL;
  ObitFArray   *ConvFn=NULL;
  ObitImage    *scrImage=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM], blc0[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM], trc0[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar        *today=NULL;
  olong        i, pos[IM_MAXDIM], Cen[2], temp, naxis[2];
  gboolean     odd;
  gchar        *routine = "ObitImageUtilVel";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));

  /* Control parameters */
  ObitInfoListGetP(inImage->info, "Parms", &type, dim, (gpointer)&Parms); 
  for (i=0; i<IM_MAXDIM; i++) blc[i] = blc0[i];
  ObitInfoListGetTest(inImage->info, "BLC", &type, dim, blc);
  for (i=0; i<IM_MAXDIM; i++) trc[i] = trc0[i];
  ObitInfoListGetTest(inImage->info, "TRC", &type, dim, trc);

  /* Open and close to get full size */
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (inImage->info, "BLC", OBIT_long, dim, blc0);
  ObitInfoListAlwaysPut (inImage->info, "TRC", OBIT_long, dim, trc0);
  iretCode = ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  iretCode = ObitImageClose (inImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Set actual requested BLC, TRC */
  ObitInfoListAlwaysPut (inImage->info, "BLC", OBIT_long, dim, blc);
  for (i=0; i<3; i++) if (trc[i]<=0.0) trc[i] = inImage->myDesc->inaxes[i];
  /* First two must be even for convolution to work */
  for (i=0; i<2; i++) {
    temp = trc[i] - blc[i] + 1;
    odd = (2*(temp/2)) != temp;
    if (odd) { /* make size even - fiddle TRC is needed */
      if (trc[i]<inImage->myDesc->inaxes[i]) trc[i] += 1.0;
      else trc[i] -= 1.0;
    }
  }
  ObitInfoListAlwaysPut (inImage->info, "TRC", OBIT_long, dim, trc);

  /* Blank or zero for out of range points? */
  if (Parms[2]<=0.5) newVal = fblank;
  else newVal = 0.0;

  /* Make scratch copy for convolution */
  scrImage = newObitImageScratch (inImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Open input image */
  iretCode = ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Copy descriptor */
  outImage->myDesc = ObitImageDescCopy(inImage->myDesc, outImage->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Setup for convolutions - make convolving function */
  naxis[0] = ObitFFTSuggestSize(inImage->myDesc->inaxes[0]);
  naxis[1] = ObitFFTSuggestSize(inImage->myDesc->inaxes[1]);
  ConvFn = ObitFArrayCreate("ConvFn", 2, naxis);
  if (Parms[3]<=0.0) FWHM = 5.0;
  else FWHM   = Parms[3];
  Cen[0] = naxis[0]/2; 
  Cen[1] = naxis[1]/2; 
  ObitFArray2DCGauss (ConvFn, Cen, FWHM);

  /* Make accumulation arrays */
  accum = newObitFArray("Accumulator");
  ObitFArrayClone (inImage->image, accum, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  ObitFArrayFill (accum, fblank);
  accumWt = newObitFArray("Accumulator");
  ObitFArrayClone (inImage->image, accumWt, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  ObitFArrayFill (accumWt, fblank);

  /* Close input for convolution */
  iretCode = ObitImageClose (inImage, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  
  /* Convolve to scratch image */
  ObitConvUtilConv (inImage, ConvFn, FALSE, 1.0, scrImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  
  /* Creation date today */
  today = ObitToday();
  strncpy (outImage->myDesc->date, today, IMLEN_VALUE-1);
  if (today) g_free(today);
 
  /* But only one plane */
  outImage->myDesc->inaxes[2] = 1;

  /* Force to float pixels */
  outImage->myDesc->bitpix=-32;

   /* (re)Open input image */
  iretCode = ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  
 /* Open scratch/convolved image */
  iretCode = ObitImageOpen (scrImage, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, scrImage->name);

  /* Loop over planes until hitting EOF */
  iplane = 0;
  while (iretCode!= OBIT_IO_EOF) {
    iplane++;
    /* Read Image */
    iretCode = ObitImageRead (inImage, NULL, err);
    if (iretCode == OBIT_IO_EOF) break;
    if (err->error) Obit_traceback_msg (err, routine, inImage->name);

    /* Read Convolved Image */
    iretCode = ObitImageRead (scrImage, NULL, err);
    if (iretCode == OBIT_IO_EOF) break;
    if (err->error) Obit_traceback_msg (err, routine, scrImage->name);

    /* Get plane Max value */
    maxF = fabs (ObitFArrayMaxAbs(inImage->image, pos));
    
    /* RMS of convolved image */
    RMS  = ObitFArrayRMS(scrImage->image);
			 
    /* Blank out of range pixels in convolved image */
    minAllow = MAX (Parms[0]*RMS, Parms[1]*maxF);
    ObitFArrayInClip (scrImage->image, -minAllow, minAllow, fblank);

    /* Blank Image plane */
    ObitFArrayBlank (inImage->image, scrImage->image, inImage->image);

    /* Accumulate Flux as Weights */
    ObitFArraySumArr (inImage->image, accumWt, accumWt);

    /* Channel Velocity */
    vel = inImage->myDesc->crval[2] + inImage->myDesc->cdelt[2] * 
      (iplane - inImage->myDesc->crpix[2]);
    
    /* Multiply image by Velocity */
    ObitFArraySMul (inImage->image, vel);
    
    /* Accumulate */
    ObitFArraySumArr (inImage->image, accum, accum);
    
  } /* End loop over input image */

  /* Close input */
  iretCode = ObitImageClose (inImage, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  /* Free image buffer */
  inImage->image = ObitFArrayUnref(inImage->image);

  /* Close scratch */
  iretCode = ObitImageClose (scrImage, err);
  if (err->error) Obit_traceback_msg (err, routine, scrImage->name);
  scrImage = ObitImageUnref(scrImage); /* release scratch image */

  /* Normalize image */
  ObitFArrayDiv (accum, accumWt, accum);

  /* newVal Replace blanks if needed */
  if (newVal != fblank)
    ObitFArrayDeblank (accum, newVal);


  /* Open output image */
  /* Use external buffer for writing output */
  outImage->extBuffer = TRUE;
  oretCode = ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  /* Write plane */
  oretCode = ObitImageWrite(outImage, accum->array, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  oretCode = ObitImageClose (outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  /* Unset external buffer for writing */
  outImage->extBuffer = FALSE;
  accum   = ObitFArrayUnref(accum);
  accumWt = ObitFArrayUnref(accumWt);
  ConvFn  = ObitFArrayUnref(ConvFn);

} /* end  ObitImageUtilVel */

/**
 * Copy an image with selection by BLC, TRC, inc
 * \param inImage  Input image
 * Parameters in info:
 * \li "BLC"     OBIT_int (7) Lowest x,y,v pixel number selected [def 1,1,1]
 * \li "TRC"     OBIT_int (7) Highest x,y,v pixel number selected [def all]
 * \li "inc"     OBIT_int (7) Pixel increment on each axis [def all 1]
 * \param outImage Image to be written.  Must be previously instantiated.
 * \param err      Error stack, returns if not empty.
 */
void 
ObitImageUtilSelCopy (ObitImage *inImage, ObitImage *outImage, ObitErr *err)
{
  ObitIOCode   iretCode, oretCode;
  ofloat       tmp;
  olong         iplane,itemp,  plane[IM_MAXDIM-2] = {0,1,1,1,1};
  olong        lblc[IM_MAXDIM], ltrc[IM_MAXDIM], linc[IM_MAXDIM];
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         inc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         blc[IM_MAXDIM], blc0[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM], trc0[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar        *today=NULL;
  olong        i;
  gboolean     want;
  gchar        *routine = "ObitImageUtilSelCopy";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));

  /* Control parameters */
  ObitInfoListGetTest(inImage->info, "inc", &type, dim, inc); 
  for (i=0; i<IM_MAXDIM; i++) {
    if (inc[i]<=0) inc[i] = 1;
    linc[i] = inc[i];
  }
  for (i=0; i<IM_MAXDIM; i++) blc[i] = blc0[i];
  ObitInfoListGetTest(inImage->info, "BLC", &type, dim, blc);
  for (i=0; i<IM_MAXDIM; i++) trc[i] = trc0[i];
  ObitInfoListGetTest(inImage->info, "TRC", &type, dim, trc);

  /* Open and close to get full size */
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (inImage->info, "BLC", OBIT_long, dim, blc0);
  ObitInfoListAlwaysPut (inImage->info, "TRC", OBIT_long, dim, trc0);
  iretCode = ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  iretCode = ObitImageClose (inImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Set actual requested BLC, TRC */
  ObitInfoListAlwaysPut (inImage->info, "BLC", OBIT_long, dim, blc);
  ObitInfoListAlwaysPut (inImage->info, "TRC", OBIT_long, dim, trc);

  /* Open input image */
  iretCode = ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Copy descriptor */
  outImage->myDesc = ObitImageDescCopy(inImage->myDesc, outImage->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Creation date today */
  today = ObitToday();
  strncpy (outImage->myDesc->date, today, IMLEN_VALUE-1);
  if (today) g_free(today);
 
  /* Output image size etc. */
  for (i=0; i<IM_MAXDIM; i++) {
    tmp = 0.99 + (ofloat)inImage->myDesc->inaxes[i] / (ofloat)inc[i];
    outImage->myDesc->inaxes[i] = (olong)tmp;
    outImage->myDesc->cdelt[i]  = inImage->myDesc->cdelt[i] * inc[i];
    outImage->myDesc->crpix[i]  = 1.0 + 
      ((inImage->myDesc->crpix[i]-1.0) / inc[i]);
    /* Reading will subimage array set blc, trc to whole array */
    lblc[i] = 1;
    ltrc[i] = inImage->myDesc->inaxes[i];
  }

  /* Force to float pixels */
  outImage->myDesc->bitpix=-32;

  /* Open/create output image */
  oretCode = ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Loop over planes until hitting EOF */
  iplane = 0;
  while (iretCode!= OBIT_IO_EOF) {
    /* Which plane */
    plane[0]++;
    for (i=0; i<IM_MAXDIM-3; i++) {
      if (plane[i] > inImage->myDesc->inaxes[i-2]) {
	plane[i+1]++;
	plane[i] = 1;
      } else break;
    } /* end loop updating planes */

    iretCode = ObitImageRead (inImage, NULL, err);
    if (iretCode == OBIT_IO_EOF) break;
    if (err->error) Obit_traceback_msg (err, routine, inImage->name);

    /* This plane selected? */
    want = TRUE;
    for (i=0; i<IM_MAXDIM-3; i++) {
      itemp = 1 + (plane[i]-1) / inc[i+2];
      want = want && (itemp*inc[i+2] == plane[i]);
    }
    if (!want) continue;

    /* Select pixels on plane */
    ObitFArraySelInc (inImage->image, outImage->image, lblc, ltrc, linc, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);

    /* Write plane */
    oretCode = ObitImageWrite(outImage, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    
  } /* End loop over input image */

  /* Close files */
  iretCode = ObitImageClose (inImage, err);
  oretCode = ObitImageClose (outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  /* Free image buffers */
  inImage->image  = ObitFArrayUnref(inImage->image);
  outImage->image = ObitFArrayUnref(outImage->image);

} /* end  ObitImageUtilSelCopy */

/**
 * Filter an image outside of a radius from the origin in FT space.
 * Intended to filter out out of band noise in single dish images.
 * Filters by a function with 1.0/(nx*ny) inside radius and outside tapers
 * by an exponential with scale distance 10 pixels.
 * \param inImage  Input Image
 * \param outImage Output image, may be inImage
 * \param radius   distance from origin in uv space (m)
 * \param err      Error stack, returns if not empty.
 */
void ObitImageUtilUVFilter (ObitImage *inImage, ObitImage *outImage, ofloat radius, 
			    ObitErr *err)
{
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOSize IOBy;
  gchar *today=NULL;
  ofloat Lambda, dx, dy, dist, xcenter, ycenter, val, radius2, arg;
  olong FFTdim[2], cen[2], ix, iy;
  ObitFArray *inFArray=NULL, *outFArray=NULL, *inFArrayCopy=NULL, *maskArray=NULL;
  ObitCArray *inCArray=NULL, *outCArray=NULL;
  ObitFFT *FFTfor=NULL, *FFTrev=NULL;
  gchar *routine = "ObitImageUtilUVFilter";
    
  /* error checks */
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));

 /* Do I/O by plane and all of plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListAlwaysPut (inImage->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy);
  ObitInfoListAlwaysPut (outImage->info,"IOBy", OBIT_long, dim, (gpointer)&IOBy);
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (inImage->info, "BLC",  OBIT_long, dim, blc); 
  ObitInfoListAlwaysPut (inImage->info, "TRC",  OBIT_long, dim, trc);

  /* Open input */
  ObitImageOpen(inImage, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Read input plane */
  ObitImageRead (inImage, NULL, err) ;
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* FFT size */
  FFTdim[0] = ObitFFTSuggestSize(inImage->myDesc->inaxes[0]);
  FFTdim[1] = ObitFFTSuggestSize(inImage->myDesc->inaxes[1]);

  /* Create float arrays for FFT size */
  inFArray  = ObitFArrayCreate("input",  2, FFTdim);
  outFArray = ObitFArrayCreate("output", 2, FFTdim);
    
  /* Save input for blanking in case overwritten */
  if (ObitImageSame(inImage, outImage, err)) {
    inFArrayCopy =  ObitFArrayCopy(inImage->image, inFArrayCopy, err);
  } else inFArrayCopy =  ObitFArrayRef(inImage->image);

  /* Pad input into work FArray */
  ObitFArrayPad(inImage->image, inFArray, 1.0);
  /* and God said "The center of an FFT will be at the corners"*/
  ObitFArray2DCenter (inFArray);
  /* Zero output FArray and use as imaginary part */
  ObitFArrayFill(outFArray, 0.0);
  
  /* Create FFT for full complex FFT */
  FFTfor = newObitFFT("Forward FFT", OBIT_FFT_Forward, OBIT_FFT_FullComplex, 2, FFTdim);
    
  /* Create complex arrays for FFT size */
  inCArray  = ObitCArrayCreate("input", 2, FFTdim);
  outCArray = ObitCArrayCreate("output", 2, FFTdim);
    
  /* Copy input to scratch CArray */
  ObitCArrayComplex(inFArray, outFArray, inCArray);

  /* close input */  
  ObitImageClose (inImage, err) ;
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Forward FFT */
  ObitFFTC2C (FFTfor, inCArray, outCArray);
    
  /* Create aperature mask */
  maskArray = ObitCArrayMakeF (outCArray);
  /* Scaling for FFT */
  ObitFArrayFill (maskArray, 1.0/(FFTdim[0]*FFTdim[1]));
  /* Wavelength */
  Lambda = 2.997924562e8 / inImage->myDesc->crval[inImage->myDesc->jlocf];
  /* Pixel size in uv plane in m */
  dx = fabs((FFTdim[0]*0.25*inImage->myDesc->cdelt[0]/57.296) / Lambda);
  dy = fabs((FFTdim[1]*0.25*inImage->myDesc->cdelt[1]/57.296) / Lambda);

  /* Form mask */
  xcenter = (ofloat)(FFTdim[0]/2); 
  ycenter = (ofloat)(FFTdim[1]/2);
  radius2 = radius*radius;
  for (iy=0; iy<FFTdim[1]; iy++) {
    for (ix=0; ix<FFTdim[1]; ix++) {
      dist = dx*(ix-xcenter)*dx*(ix-xcenter) + dy*(iy-ycenter)*dy*(iy-ycenter);
      /* Add exp taper outside of radius */
      if (dist>radius2) {
	dist = sqrt(dist);
	val = maskArray->array[iy*FFTdim[0]+ix];
	arg = -(dist-radius)/(10.0*dx);
	if (arg>-14.0) val *= exp(arg);
	else val = 0.0;
	maskArray->array[iy*FFTdim[0]+ix] = val;
      } /* End outside radius */
    } /* end loop in x */
  } /* end loop in y */
  
  /* Mask*/
  ObitFArray2DCenter(maskArray); /* Swaparonie */
  ObitCArrayFMul (outCArray, maskArray, outCArray);

  /* Create FFT for back FFT */
  FFTrev = newObitFFT("Forward FFT", OBIT_FFT_Reverse, OBIT_FFT_FullComplex, 2, FFTdim);

  /* Back FFT */
  ObitFFTC2C(FFTrev, outCArray, inCArray);
    
  /* Extract Real */
  ObitCArrayReal (inCArray, outFArray);
  /* and God said "The center of an FFT will be at the corners" */
  ObitFArray2DCenter(outFArray);
    
  /* Setup output */
  /* If input not output, clone */
  if (!ObitImageSame(inImage, outImage, err)) {
    ObitImageClone (inImage, outImage, err);
  }
  ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  /* Creation date today */
  today = ObitToday();
  strncpy (outImage->myDesc->date, today, IMLEN_VALUE-1);
  if (today) g_free(today);
  
  /* Extract output portion */
  cen[0] = FFTdim[0]/2; cen[1] = FFTdim[1]/2;
  blc[0] = cen[0] - inFArrayCopy->naxis[0] / 2; 
  trc[0] = cen[0] - 1 + inFArrayCopy->naxis[0] / 2;
  blc[1] = cen[1] - inFArrayCopy->naxis[1] / 2; 
  trc[1] = cen[1] - 1 + inFArrayCopy->naxis[1] / 2;
  ObitImageUnref(outImage->image);
  outImage->image = ObitFArraySubArr(outFArray, blc, trc, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Blank output where input blanked */
  ObitFArrayBlank (outImage->image, inFArrayCopy, outImage->image);

  /* Write */
  ObitImageWrite (outImage, outImage->image->array, err);
  ObitImageClose (outImage,err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Cleanup */
  if (inImage->image)  inImage->image  = ObitFArrayUnref(inImage->image);
  if (outImage->image) outImage->image = ObitFArrayUnref(outImage->image);
  if (inFArray)     ObitFArrayUnref(inFArray);
  if (outFArray)    ObitFArrayUnref(outFArray);
  if (inFArrayCopy) ObitFArrayUnref(inFArrayCopy);
  if (maskArray)    ObitFArrayUnref(maskArray);
  if (inCArray)     ObitCArrayUnref(inCArray);
  if (outCArray)    ObitCArrayUnref(outCArray);
  if (FFTfor)       ObitFFTUnref(FFTfor);
  if (FFTrev)       ObitFFTUnref(FFTrev);
} /* end ObitImageUtilUVFilter */

/*----------------------Private functions---------------------------*/

/**
 * Fills an existing character array with the string for the current date.
 * \param date Character string to accept the string (10 char+null)
 * \param len  Actual length of date (should be at least 11)
 */
static void ObitImageUtilCurDate (gchar *date, olong len)
{
  struct tm *lp;
  time_t clock;
  
  /* Get time since 00:00:00 GMT, Jan. 1, 1970 in seconds. */
  time (&clock);
  
  /* Convert to  broken-down time. */
  lp = localtime (&clock);
  lp->tm_mon++; /* For some bizzare reason, month is 0-rel */
  
  /* Full year */
  if (lp->tm_year<1000)  lp->tm_year += 1900; 
  lp->tm_mon = MAX (1, lp->tm_mon);

  /* to output */
  g_snprintf (date, len, "%4.4d-%2.2d-%2.2d",
	      lp->tm_year, lp->tm_mon, lp->tm_mday);
} /* end ObitImageUtilCurDate */

/**
 * Fit Gaussian to a Beam.
 * Adapted from AIPS BMSHP.FOR, FITBM.FOR
 * \param beam  Beam image to fit.
 * \param err   Error stack, returns if not empty.
 */
static void ObitImageUtilFitBeam (ObitImage *beam, ObitErr *err)
{
  olong blc[2], trc[2], center[2], i, j, k, l, ijk, iflip, irow, ilast, naxis[2];
  olong ierr;
  ofloat peak, *array, dx, dy, *x;
  ofloat  p[3], y[3];
  ofloat cellx, celly, xfact, c1, c2, c3, temp, s, c, ss, cc;
  ofloat bmaj, bmin, bpa;
  ObitFArray *beamData = NULL, *beamCenter = NULL, *matx = NULL;
  gchar *routine = "ObitImageUtilFitBeam";

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(beam));

  /* Trim edges of beam */
  blc[0] = 5;
  blc[1] = 5;
  trc[0] = beam->image->naxis[0]-5;
  trc[1] = beam->image->naxis[1]-5;
  beamCenter = ObitFArraySubArr (beam->image, blc, trc, err);
  if (err->error) goto cleanup;

  /* Find center */
  peak = ObitFArrayMax (beamCenter, center);

  /* Center better be 1.0 */
  if (fabs (peak-1.0) > 0.001) { 
    /* Close so the image can be looked at */
    if ((ObitImageClose (beam, err)  
	 != OBIT_IO_OK) || (err->error>0)) {  /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s",  
      		     routine, beam->name); 
      goto cleanup; 
    }
    Obit_log_error(err, OBIT_Error,  
      		   "Beam peak (%f) not 1.0 for %s",  
      		   peak, beam->name); 
    goto cleanup; 
  } 

  /* extract half center into an ObitFArray */
  blc[0] = center[0] - 5;
  blc[1] = center[1];
  trc[0] = center[0] + 5;
  trc[1] = center[1] + 5;
  beamData = ObitFArraySubArr (beamCenter, blc, trc, err);
  if (err->error) goto cleanup;

  /* Image info */
  cellx = beam->myDesc->cdelt[0];
  celly = beam->myDesc->cdelt[1];
  xfact = fabs (cellx);

  /* Make an array for the solution matrix */
  naxis[0] = 3; naxis[1] = 3;
  matx = ObitFArrayCreate ("Soln Matrix", 2, naxis);
  x = matx->array;
  for (i=0; i<9; i++) x[i] = 0.0;

  /* Init */
  y[0] = y[1] = y[2] = 0.0;

  /* Loop down rows  doing alternate halves. go only to first
     decending 0.35 from center. */
  iflip = 1;
  array = beamData->array;
  for (i=0; i<6; i++) {  /* loop over rows */
    for (ijk= 0; ijk<2; ijk++) { /* loop 65 */
      iflip = -iflip;
      ilast =  5 - iflip;
      
      /* Loop over row flipping halves */
      for (j= ijk; j<6; j++) { /* loop 60 */
	irow = 5 + j * iflip; /* which cell? */

	/* Have we gone far enough? */
	if ((array[irow] < 0.35)  &&  
	    (array[irow] <  array[ilast])) break;
	if (array[irow] < 0.35) continue;
	ilast = irow;

	/* compute displacements from center */
	dx = iflip * j * cellx / xfact;
	dy = -i * celly / xfact;

	/* compute partials wrt c1,c2,c3 */
	p[0] = dx * dx;
	p[1] = dy * dy;
	p[2] = dx * dy;

	/* sum partials into x matrix and y vector */
	for (k= 0; k<3; k++) { 
	  y[k] = y[k] - log (array[irow]) * p[k];
	  for (l= 0; l<3; l++) { 
	    x[l*3 + k] += p[k] * p[l];
	  } 
	}
      } /* end loop  L60:  */;
    } /* end loop  L65:  */;
    
    array += beamData->naxis[0]; /* next row in array */
  }/* end loop over rows */

  /* Do beam fitting use ObitFArray to invert matrix */
  ObitFArray2DSymInv(matx, &ierr);
  if (ierr!=0) { /* fit failed issue warning and use 1 cell beam */
   Obit_log_error(err, OBIT_InfoWarn, 
		  "Solution for Restoring beam failed for %s", 
		  beam->name);
   /* default solution */
   c1 = 0.5;
   c2 = 0.5;
   c3 = 0.0;
  } else { 
    /* actual solution */
    c1 = y[0] * x[0*3+0] + y[1] * x[0*3+1] + y[2] * x[0*3+2];
    c2 = y[0] * x[1*3+0] + y[1] * x[1*3+1] + y[2] * x[1*3+2];
    c3 = y[0] * x[2*3+0] + y[1] * x[2*3+1] + y[2] * x[2*3+2];
  }

  /* convert to sigmas and pa. */
  /* make sure arg of atan is determinate. */
  if (fabs (c1-c2) < 1.0e-10) {
    bpa = 45.0;
    if (c3<0) bpa = -bpa;
  } else {
    bpa = 28.6478 * atan (c3 / (c1-c2));
  } 
  
  /* compute sigma**2 */
  s = sin (bpa/57.29578);
  c = cos (bpa/57.29578);
  ss = s * s;
  cc = c * c;
  if (fabs (ss-cc)  >  1.0e-4) {
    bmaj = 0.5 * (cc*cc - ss*ss) / (cc*c1 - ss*c2);
    bmin = 0.5 * (ss*ss - cc*cc) / (ss*c1 - cc*c2);
    /* special case near 45 deg. */
  } else {
    bmaj = c3 / (4.0 * s * c) + c1;
    bmin = 0.5 / (c1/cc - bmaj);
    bmaj = 0.5 / bmaj;
  } 

  /* check if soln. is real */
  /* use default */
  if ((bmaj <= 0.0)  ||  (bmin <= 0.0)) {
    bmaj = 1.0;
    bmin = 1.0;
    bpa = 0.0;
    Obit_log_error(err, OBIT_InfoWarn, 
		   "Solution for Restoring beam failed for %s", 
		   beam->name);
  } 

  /* convert to sigmas. */
  bmaj = sqrt (bmaj) * xfact;
  bmin = sqrt (bmin) * xfact;
  /* up to here bmaj is minor axis: rest of pgm wants it as  major 
     - fix this here */
  if (bmaj > bmin) {
    bpa = bpa - 90.0;
  } else {
    temp = bmaj;
    bmaj = bmin;
    bmin = temp;
  } 

  /* convert to fwhm */
  bmaj = bmaj * 2.3548;
  bmin = bmin * 2.3548;

  /* add map rotation, force to +/- 90 deg. */
  bpa = bpa - beam->myDesc->crota[1];
  if (bpa > 90.0) bpa = bpa - 180.0;
  if (bpa < -90.0) bpa = bpa + 180.0;
  
  /* Give informative message about fitted size */
  Obit_log_error(err, OBIT_InfoErr, 
		 "Fitted beam %f x %f asec, PA = %f for %s", 
		 bmaj*3600.0, bmin*3600.0, bpa, beam->name);

  /* Save fitted values */
  beam->myDesc->beamMaj = bmaj;
  beam->myDesc->beamMin = bmin;
  beam->myDesc->beamPA  = bpa;

  /* Cleanup */
 cleanup: matx = ObitFArrayUnref (matx);
  beamData   = ObitFArrayUnref(beamData);
  beamCenter = ObitFArrayUnref(beamCenter);
 
  /* Free image buffer */
  beam->image = ObitFArrayUnref(beam->image);
  if (err->error) Obit_traceback_msg (err, routine, beam->name);

} /* end ObitImageUtilFitBeam */

