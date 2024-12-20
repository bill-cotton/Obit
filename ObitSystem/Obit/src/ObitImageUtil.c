/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2023                                          */
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
#include "ObitGPUGrid.h"
#include "ObitUVGrid.h"
#include "ObitImageDesc.h"
#include "ObitThread.h"
#include "ObitUVGridWB.h"
#include "ObitUVGridMF.h"
#include "ObitImageWB.h"
#include "ObitImageMF.h"
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
#include "ObitUVImager.h"
#include "ObitFFT.h"
#include "ObitTableCCUtil.h"
#include "ObitImageFit.h"
#include "ObitBeamShape.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageUtil.c
 * ObitImageUtil module function definitions for ObitImage class.
 */

/*---------------Private structures----------------*/
/* Image interpolation threaded function argument */
typedef struct {
  /* Input descriptor */
  ObitImageDesc *inDesc;
  /* Input plane pixel data */
  ObitFArray *inData;
  /* Output descriptor */
  ObitImageDesc *outDesc;
  /* Output plane pixel data */
  ObitFArray *outData;
  /* X Pixel array */
  ObitFArray *XPixData;
  /* Y Pixel array */
  ObitFArray *YPixData;
  /* Also do Weights? */
  gboolean   doWeight;
  /* Output weight plane pixel data */
  ObitFArray *wtData;
  /* Radius in pixels of weighting circle */
  olong      radius;
  /* Number of Zernike corrections */
  olong      nZern;
  /* Zernike coefficients */
  ofloat     *ZCoef;
  /* First (1-rel) row in image to process this thread */
  olong      first;
  /* Highest (1-rel) row in image to process this thread  */
  olong      last;
  /* thread number, <0 -> no threading  */
  olong      ithread;
  /* Obit Thread object */
  ObitThread  *thread;
  /* Obit error stack object */
  ObitErr    *err;
  /* Input Image Interpolator */
  ObitFInterpolate *Interp;
} InterpFuncArg;

/* Image primary beam  threaded function Argument */
typedef struct {
  /** First (1-rel) value in y to process this thread */
  olong         first;
  /** Highest (1-rel) value in y to process this thread  */
  olong         last;
  /** thread number, >0 -> no threading   */
  olong         ithread;
  /** Output descriptor */
  ObitImageDesc *outDesc;
  /** Output beam image data */
  ObitFArray    *outArr;
  /** Beam Shape */
  ObitBeamShape *bs;
  /** Invert? */
  gboolean      doInvert;
  /** Obit Thread object */
  ObitThread  *thread;
  /* Obit error stack object */
  ObitErr    *err;
} PBCorFuncArg;
/*---------------Private function prototypes----------------*/
/** Private: Get Date string for current date */
static void ObitImageUtilCurDate (gchar *date, olong len);

/** Private: Threaded Image interpolator */
static gpointer ThreadImageInterp (gpointer arg);

/** Private: Threaded Calculate pixels */
static gpointer ThreadGetXYPixels (gpointer arg);

/** Private: Make Threaded Image interpolator args */
static olong MakeInterpFuncArgs (ObitThread *thread, olong radius,
				 olong nZern, ofloat *ZCoef,
				 ObitImageDesc *inDesc, ObitImageDesc *outDesc, 
				 ObitFInterpolate *Interp,
				 ObitErr *err, InterpFuncArg ***ThreadArgs);

/** Private: Delete Threaded Image interpolator args */
static void KillInterpFuncArgs (olong nargs, InterpFuncArg **ThreadArgs);

/** Private: Threaded Image interpolator */
static gpointer ThreadImagePBCor (gpointer arg);

/** Private: Make Threaded Image interpolator args */
static olong MakePBCorFuncArgs (ObitThread *thread,
				ObitBeamShape *bs, ObitImageDesc *outDesc, 
				ObitFArray *outArr, gboolean doInvert,
				ObitErr *err, PBCorFuncArg ***ThreadArgs);

/** Private: Delete Threaded Image PBCor args */
static void KillPBCorFuncArgs (olong nargs, PBCorFuncArg **ThreadArgs);

/** Private: Set convolution kernal for 3rd axis */
static void SetConvKernal3 (ofloat Target, olong npln, olong hwidth, 
			    olong *Start, ofloat *Kernal);


/*----------------------Public functions---------------------------*/

/**
 * Create basic ObitImage structure and fill out descriptor.
 * Imaging parameters are on the inUV info member as arrays for a number 
 * of fields.
 * \li "nChAvg" OBIT_long (1,1,1) number of channels to average.
 *              This is for spectral line observations and is ignored
 *              if the IF axis on the uv data has more than one IF.
 *              Default is continuum = average all freq/IFs. 0=> all.
 * \li "rotate" OBIT_float (?,1,1) Desired rotation on sky (from N thru E) in deg. [0]
 * \li "nx"     OBIT_long (?,1,1) Dimension of image in RA [no default].
 *              This and the following are arrays with one entry per field.
 * \li "nxBeam" OBIT_long (?,1,1) Dimension of beam in RA, [def. nx]
 * \li "ny"     OBIT_long (?,1,1) Dimension of image in declination[no default]
 * \li "nyBeam" OBIT_long (?,1,1) Dimension of beam in declination, [def. ny]
 * \li "xCells" OBIT_float (?,1,1) X (=RA) cell spacing in asec [no default]
 * \li "yCells" OBIT_float (?,1,1) Y (=dec) cell spacing in asec [no default]
 * \li "xShift" OBIT_float (?,1,1) Desired shift in X (=RA) in degrees. [0]
 * \li "yShift" OBIT_float (?,1,1) Desired shift in Y (=dec) in degrees. [0]
 * \li "do3D"   OBIT_bool (1,1,1) 3D image, else 2D? [def TRUE]
 * \li "doGrid" OBIT_bool (?,1,1) if do3D=FALSE, force image to common grid [TRUE]
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
  gboolean *barray, doGrid;
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

  /* Align 2D to grid? */
  doGrid = TRUE;
  ObitInfoListGetP(inUV->info, "doGrid", &type, dim, (gpointer*)&barray);
  if (barray!=NULL) doGrid = barray[fieldNo-1];

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

  /* 2D/3D imaging */
  outImage->myDesc->do3D = TRUE;
  ObitInfoListGetTest(inUV->info, "do3D", &type, dim, &outImage->myDesc->do3D);

  /* Fill in descriptor */
  ObitImageUtilUV2ImageDesc (inUV->myDesc, outImage->myDesc, doGrid, nChAvg);

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
   theBeam->myDesc->do3D      = outImage->myDesc->do3D;
   /* Fill in descriptor */
   ObitImageUtilUV2ImageDesc (inUV->myDesc, theBeam->myDesc, FALSE, nChAvg);
   /* Header may have changed */
   if (theBeam->myStatus!=OBIT_Inactive) theBeam->myStatus = OBIT_Modified;
 }

  /* Replace uv data shift */
  inUV->myDesc->xshift = xShift;
  inUV->myDesc->yshift = yShift;

  return outImage;
} /* end ObitImageUtilCreateImage */

/**
 *  Make an image from from uv data, info in ObitInfoList.
 * Grids, FFTs and makes corrections for the gridding convolution.
 * This interface allows multi-processing and/or multi-threading.
 * \param inList  Input File InfoList
 *  Input UV data prefix = "ImgUV"
 *      Input uv data. Should be in form of stokes to be imaged
 *      will all calibration and selection applied and any 
 *      weighting applied.
 * Weighting/Imaging parameters on inList copied to uvdata:
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
 * \li "do3D"   OBIT_bool (1,1,1) 3D image, else 2D? [def TRUE]
 * \li "WtPower" OBIT_float scalar = Power to raise weights to.  [def = 1.0]
 *              Note: a power of 0.0 sets all the output weights to 1 as modified
 *              by uniform/Tapering weighting.  Applied in determining weights 
 *              as well as after.
 * 
 *  Output Image prefix = "OutImg" Image to be written.  Must be previously defined.
 *      Beam normalization factor is written to output Beam
 *      infoList as SUMWTS
 *    "doBeam"  gboolean if TRUE also make beam.  Will make the myBeam member of [Def F]
 *               outImage.
 *               If FALSE, and myGrid->BeamNorm 0.0 then reads SUMWTS value 
 *               from beam infolist
 *    "doWeight" gboolean if TRUE Apply uniform weighting corrections to uvdata
 *                 before imaging. [Def F]
 *    "doFlatten" gboolean if TRUE Flatten mosaic when done [Def F]
 *    "field" olong  Which field (1-rel) to Image, 0=> all [def 0]
 *
 *  Output Beam prefix (if doBeam) = "OutBeam" Image to be written.  
 *      Must be previously defined.
 *      If doBeam=FALSE, the following must be given in inList:
 *     "SUMWTSnnnn" float Beam normalization factor where nn is the 1-rel field number - 
 *      returned if doBeam=TRUE
 * \param err      Error stack, returns if not empty.
 */
void ObitImageUtilMakeImageFileInfo (ObitInfoList *inList, ObitErr *err)
{
  ObitUVImager *myUVImager = NULL;
  gboolean doBeam, doWeight, doFlatten;
  ObitImage *theBeam;
  olong field[1], ifield, lofield, hifield;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar strsumwt[16];
  ofloat sumwts;
  ObitUVImagerClassInfo *imgClass;
  gchar *routine = "ObitImageUtilMakeImageFileInfo";

  /* error checks */
  if (err->error) return;

  /* Control */
  doBeam  = FALSE;
  ObitInfoListGetTest (inList, "doBeam",   &type, dim, &doBeam); 
  doWeight  = FALSE;
  ObitInfoListGetTest (inList, "doWeight", &type, dim, &doWeight); 
  doFlatten = FALSE;
  ObitInfoListGetTest (inList, "doFlatten", &type, dim, &doFlatten); 
  field[0]     = 0;
  ObitInfoListGetTest (inList, "field", &type, dim, field); 

  /* Build object */
  myUVImager = ObitUVImagerFromInfo("ImgUV", inList, err);
  if (err->error) Obit_traceback_msg (err, routine, myUVImager->name);

  /* Which fields */
  if (field>0) {
    lofield = field[0];
    hifield = field[0];
  } else {
    lofield = 1;
    hifield = myUVImager->mosaic->numberImages;
  }

  /* If not making beam copy SUMWTS */
  if (!doBeam) {
    for (ifield=lofield; ifield<=hifield; ifield++) {
      sprintf (strsumwt,"SUMWTS%5.5d",ifield);
      sumwts = 0.0;
      ObitInfoListGetTest(inList, strsumwt, &type, dim, &sumwts);
      dim[0] = 1;
      if (myUVImager->mosaic->images[ifield-1]->myBeam) {
	theBeam = (ObitImage*)myUVImager->mosaic->images[ifield-1]->myBeam;
	ObitInfoListAlwaysPut(theBeam->info, "SUMWTS", OBIT_float, dim, &sumwts);
      }
    }
  }

  /* Make image */
  imgClass = (ObitUVImagerClassInfo*)myUVImager->ClassInfo;  
  imgClass->ObitUVImagerImage (myUVImager, field, doWeight, doBeam, doFlatten, err);
  if (err->error) Obit_traceback_msg (err, routine, myUVImager->name);

  /* Made beam? Copy SUMWTS */
  if (doBeam) {
    for (ifield=lofield; ifield<=hifield; ifield++) {
      sprintf (strsumwt,"SUMWTS%5.5d",ifield);
      sumwts = 0.0;
      if (myUVImager->mosaic->images[ifield-1]->myBeam) {
	theBeam = (ObitImage*)myUVImager->mosaic->images[ifield-1]->myBeam;
	ObitInfoListGetTest(theBeam->info, "SUMWTS", &type, dim, &sumwts);
      }
      dim[0] = 1;
      ObitInfoListAlwaysPut(inList, strsumwt, OBIT_float, dim, &sumwts);
    }
  }

} /* end ObitImageUtilMakeImageFileInfo */

/**
 * Grids, FFTs and makes corrections for the gridding convolution.
 * Uses (creating if necessary) the myGrid member of out.
 * Can make multichannel continuum or spectral line images
 * Supports multiple GPUs
 * \param inUV     Input uv data. Should be in form of stokes to be imaged
 *                 will all calibration and selection applied and any 
 *                 weighting applied.
 * \param outImage Image to be written.  Must be previously instantiated.
 *                 Beam normalization factor is written to output Beam
 *                 infoList as SUMWTS
 *                 Supports ObitImage, bitImageMF or ObitImageWB.
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
 *              by uniform/Tapering weighting.  Applied in determing weights 
 *              as well as after.
 * \li "do3D"   OBIT_bool (1,1,1) 3D image, else 2D? [def TRUE]
 * \li "doGPUGrid" OBIT_bool (1,1,1) use GPU for Gridding [def FALSE]
 * \li "GPU_no"    OBIT_long list of GPU device numbers 0-rel [def (0)]
 * \li "chDone"    OBIT_bool (nSpec,1,1) Array of flags indicating channels are done [all FALSE]
 * \param channel  Which frequency channel to image, 0->all.
 * \param err      Error stack, returns if not empty.
 */
void ObitImageUtilMakeImage (ObitUV *inUV, ObitImage *outImage, 
			     olong channel, gboolean doBeam, gboolean doWeight,
			     ObitErr *err)
{
  ObitImage *theBeam=NULL;
  ObitImageDesc *imDesc, *bmDesc, *IODesc;
  ObitIOSize IOBy;
  ObitInfoType type;
  ofloat sumwts, imMax, imMin, BeamTaper=0., BeamNorm=0.0, Beam[3];
  gchar *outName=NULL;
  ollong *gpumem=NULL;
  olong ldevice[20], *cuda_device=NULL, iGPU=-1, num_GPU=1; 
  olong plane[5], pln, NPIO, oldNPIO, nfacet=1, ifacet=0;
  olong nplane=1, i, ichannel, icLo, icHi, norder=0, iorder=0;
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitUVGridClassInfo *gridClass;
  ObitImageClassInfo *imgClass;
  ObitUVGrid *myGrid=NULL, *beamGrid=NULL;
  gboolean doGPUGrid=FALSE, doCalSelect=FALSE, initGPU=TRUE, *chDone=NULL;
  ObitIOAccess access;
  union ObitInfoListEquiv InfoReal; 
  gchar *routine = "ObitImageUtilMakeImage";

  /* error checks */
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

  /* Get restoring beam */
  ObitInfoListGetTest(inUV->info, "Beam", &type, dim, Beam);

  /* Use GPU? */
  initGPU = TRUE;  /* Need to set up for using GPU */
  InfoReal.itg = (olong)doGPUGrid; type = OBIT_bool;
  ObitInfoListGetTest(inUV->info, "doGPUGrid", &type, (gint32*)dim, &InfoReal);
  doGPUGrid = InfoReal.itg;
  if (doGPUGrid) {
    /* How many GPUs? */
    if (ObitInfoListGetTest(inUV->info, "GPU_no", &type, (gint32*)dim, ldevice)) {
      /* Count, ignore -1s  */
      num_GPU = 0;
      for (i=0; i<dim[0]; i++) if (ldevice[i]>=0) num_GPU++;
      if (num_GPU<=0) {num_GPU=1; ldevice[0]=0;}  /* defaults to cuda device 0 */
      cuda_device = g_malloc0((num_GPU+10)*sizeof(olong));  /* Hack (+10) for sarao bug*/
      gpumem =      g_malloc0((num_GPU+10)*sizeof(ollong));
      for (i=0; i<num_GPU; i++) {
	if ((ldevice[i]<0) || (ldevice[i]>50)) ldevice[i] = i;  /* Sanity check */
	cuda_device[i] = (int)ldevice[i];
	if (ObitGPUGridCheckGPU((olong)cuda_device[i])<0) { /* Check device number */
	  Obit_log_error(err, OBIT_Error, "Invalid GPU device %d %d",cuda_device[i],i);
	  return;
	} /* end check device */
	gpumem[i] =      ObitGPUGridSetGPU (cuda_device[i]); /* initialize */
      }
    } else {  /* default GPU = 0 */
      num_GPU = 1;
      cuda_device = g_malloc0((num_GPU+10)*sizeof(olong));
      gpumem      = g_malloc0((num_GPU+10)*sizeof(ollong));
      cuda_device[0] = 0;
      if (ObitGPUGridCheckGPU((olong)cuda_device[0])<0) { /* Check device number */
	Obit_log_error(err, OBIT_Error, "Invalid GPU device %d",cuda_device[0]);
	return;
      } /* end check device */
     gpumem[0]      = ObitGPUGridSetGPU (cuda_device[0]); /* initialize */
    }
    if (err->prtLv>=2)
      Obit_log_error(err, OBIT_InfoErr, "Doing GPU Gridding with %d GPUs",num_GPU);
  } /* end using GPU */

  /* Need new gridding member? */
  IODesc = (ObitImageDesc*)outImage->myIO->myDesc;
  outName = g_strconcat ("UVGrid for: ",inUV->name,NULL);
  if (outImage->myGrid == NULL) {
    /* WB Image? */
    if (ObitImageWBIsA(outImage)) {
      outImage->myGrid = (Obit*)newObitUVGridWB(outName);
      nplane = 1;  /* Note GPU gridding doesn't work for this */
    /* MF Image? */
    } else if (ObitImageMFIsA(outImage)) {
      outImage->myGrid = (Obit*)newObitUVGridMF(outName);
      nplane = ((ObitImageMF*)outImage)->nSpec;
    /* Spectroscopic Image cube? Use ImageMF gridder */
    } else if (IODesc->inaxes[IODesc->jlocf]>1) {
      outImage->myGrid = (Obit*)newObitUVGridMF(outName);
      nplane = IODesc->inaxes[IODesc->jlocf]; 
    } else {  /* Basic image */
      outImage->myGrid = (Obit*)newObitUVGrid(outName);
      nplane = 1;
    }
    if (outName) {g_free(outName);} outName = NULL;
  }
  /* Need separate gridder for different size beam? */
  if (doBeam) {
    imDesc = outImage->myDesc;
    IODesc = ((ObitImage*)(ObitImageDesc*)outImage->myBeam)->myIO->myDesc;
    bmDesc = ((ObitImage*)outImage->myBeam)->myDesc;
    if ((imDesc->inaxes[0]!=bmDesc->inaxes[0]) || 
	(imDesc->inaxes[1]!=bmDesc->inaxes[1])) {
      outName = g_strconcat ("UVGrid for Beam: ",inUV->name,NULL);
      if (((ObitImage*)outImage->myBeam)->myGrid==NULL)
	{
	  /* WB Image? */
	  if (ObitImageWBIsA(outImage))
	    outImage->myGrid = (Obit*)newObitUVGridWB(outName);
	  /* MF Image? */
	  else if (ObitImageMFIsA(outImage))
	    outImage->myGrid = (Obit*)newObitUVGridMF(outName);
	  /* Spectroscopic Image cube? - use MF Gridding */
	  else if (IODesc->inaxes[IODesc->jlocf]>1)
	    outImage->myGrid = (Obit*)newObitUVGridMF(outName);
	  else 
	    outImage->myGrid =(Obit*) newObitUVGrid(outName);
	}
      if (outName) {g_free(outName);} outName = NULL;
    } else {  /* Same size - use image gridder */
      if (((ObitImage*)outImage->myBeam)->myGrid==NULL)
	((ObitImage*)outImage->myBeam)->myGrid = (Obit*)ObitUVGridRef(outImage->myGrid);
    }
    /* GPU Info - Beam  */
    if (doGPUGrid) {
      ifacet = 1;
      if (nfacet==1) ifacet = 0;
      iGPU = (iGPU+1)%num_GPU; /* which GPU to do this in? */
      myGrid = (ObitUVGrid*)outImage->myGrid;
      if (initGPU && (myGrid->gridInfo==NULL)) {
	myGrid->gridInfo = ObitGPUGridCreate("GPUGrid", nfacet, num_GPU, cuda_device, 
					     (Obit*)outImage, inUV, TRUE);
	initGPU = FALSE;  /* Don't do again */
      }
      iGPU = myGrid->gridInfo->cudaInfo->FacetGPU[ifacet]; /* which GPU to do this in? */
      ObitGPUGridSetGPUStruct(myGrid->gridInfo, (Obit*)myGrid, chDone, ifacet, nfacet, 
			      iGPU, num_GPU, cuda_device[iGPU], nplane, inUV, (Obit*)outImage, 
			      doBeam, err);
      if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    }
  } /* End separate beam gridder */

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
 
  /* Free Pixel array */
  outImage->image = ObitFArrayUnref(outImage->image);

  /* Get beam normalization factor from beam if needed */
  if (!doBeam) {
    /* Get Beam member from outImage */
    theBeam = (ObitImage*)outImage->myBeam;
    g_assert (ObitImageIsA(theBeam));

    sumwts = 0.0;
    ObitInfoListGetTest(theBeam->info, "SUMWTS", &type, dim, (gpointer)&sumwts);
    myGrid = (ObitUVGrid*)outImage->myGrid;
    myGrid->BeamNorm = sumwts;
    BeamNorm = sumwts;
  }    

  myGrid = (ObitUVGrid*)outImage->myGrid;
  gridClass =  (ObitUVGridClassInfo*)myGrid->ClassInfo;   /* Gridder class */
  imgClass  =  (ObitImageClassInfo*)outImage->ClassInfo;  /* Image class */

  /* Loop over channels selected */
  IODesc = (ObitImageDesc*)outImage->myIO->myDesc;
  icLo = 1; 
  icHi = IODesc->inaxes[IODesc->jlocf];
  if (channel>0) {icLo = channel; icHi = channel;}
  /* WB Image? Only 1 channel */
  if (ObitImageWBIsA(outImage)) icHi = icLo = 1;
  /* MF Image? Channels handled in gridding */
  else if (ObitImageMFIsA(outImage)) icHi = icLo = 1;
  /* Multiple parallel channels */
  else if (IODesc->inaxes[IODesc->jlocf]>1) icHi = icLo = 1;
  for (ichannel=icLo; ichannel<=icHi; ichannel++) {

    /* Make beam if needed */
    if (doBeam) {
      norder = imgClass->ObitImageGetBeamOrder(outImage); /* No. orders of beam */
      /* Loop over beams */
      for (iorder=0; iorder<=norder; iorder++) {
	/* Get relevant beam */
	theBeam = imgClass->ObitImageGetBeam(outImage, iorder, plane, err);
	if (err->error) Obit_traceback_msg (err, routine, outImage->name);
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
	
	/* Reset max/min values */
	theBeam->myDesc->maxval = -1.0e20;
	theBeam->myDesc->minval =  1.0e20;
	((ObitImageDesc*)theBeam->myIO->myDesc)->maxval = -1.0e20;
	((ObitImageDesc*)theBeam->myIO->myDesc)->minval =  1.0e20;
	
	/* Gridding setup */
	beamGrid = (ObitUVGrid*)(theBeam->myGrid);
	gridClass->ObitUVGridSetup (beamGrid, inUV, (Obit*)theBeam, TRUE, err);
	if (err->error) Obit_traceback_msg (err, routine, theBeam->name);
	
	/* Grid Beam */
	gridClass->ObitUVGridReadUV (beamGrid, inUV, err);
	if (err->error) Obit_traceback_msg (err, routine, theBeam->name);

	/* Save sum weights for higher order beams */
	if (iorder>0) beamGrid->BeamNorm = BeamNorm;
		
	/* FFT, Gridding correction, write  */
	dim[0] = dim[1] = dim[2] = 1;
	pln = ichannel;   /* Which plane to write? */
	if (ObitImageWBIsA(outImage)) pln = iorder+1;
	ObitInfoListAlwaysPut (theBeam->info, "Channel",  OBIT_long, dim,  &pln);
	gridClass->ObitUVGridFFT2Im(beamGrid, (Obit*)theBeam, err);
	if (err->error) Obit_traceback_msg (err, routine, theBeam->name);
	
	/* Use sum weights for higher order beams */
	if (iorder==0)  {
	  BeamNorm = beamGrid->BeamNorm;
	  myGrid = (ObitUVGrid*)(outImage->myGrid);
	  myGrid->BeamNorm = BeamNorm;
	}

	if (iorder==0) {  /* Stuff for dirty beam */
	  /* Tell Sum of gridding weights */
	  myGrid = (ObitUVGrid*)(outImage->myGrid);
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Sum of Weights %g for %s",myGrid->BeamNorm, outImage->name);
	  
	  /* tell Max/Min */
	  imMax = theBeam->myDesc->maxval;
	  imMin = theBeam->myDesc->minval;
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Beam max %g, min %g for %s", imMax, imMin, theBeam->name);
	  
	  /* Fit beam */
	  ObitImageUtilFitBeam (theBeam, err);
	  if (err->error) Obit_traceback_msg (err, routine, theBeam->name);

	  /* Beam specified? */
	  ObitInfoListGetTest(outImage->myDesc->info, "BeamTapr", &type, dim, &BeamTaper);
	  if ((Beam[0]>0.0) && (BeamTaper<=0.0)) {
	    Obit_log_error(err, OBIT_InfoErr, 
	     "Using Beam %f %f %f", Beam[0], Beam[1], Beam[2]);
	     theBeam->myDesc->beamMaj = Beam[0]/3600.0;
	     theBeam->myDesc->beamMin = Beam[1]/3600.0;
	     theBeam->myDesc->beamPA  = Beam[2];
	  } else if (Beam[0]<=0.0) {  /* Use first fitted beam if none given */
	    Beam[0] = theBeam->myDesc->beamMaj*3600.0;
	    Beam[1] = theBeam->myDesc->beamMin*3600.0;
	    Beam[2] = theBeam->myDesc->beamPA;
	  }
	}  /* End dirty beam stuff */
	
	/* Close Image */
	ObitImageClose (theBeam, err);
	if (err->error) Obit_traceback_msg (err, routine, theBeam->name);
	
	/* Free image buffer */
	theBeam->image = ObitFArrayUnref(theBeam->image);
      } /* end loop over beams */
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
    
    /* Reset max/min values */
    outImage->myDesc->maxval = -1.0e20;
    outImage->myDesc->minval =  1.0e20;
    ((ObitImageDesc*)outImage->myIO->myDesc)->maxval = -1.0e20;
    ((ObitImageDesc*)outImage->myIO->myDesc)->minval =  1.0e20;

    /* Gridding setup */
    gridClass->ObitUVGridSetup ((ObitUVGrid*)outImage->myGrid, inUV, (Obit*)outImage, FALSE, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    
    /* GPU Info */
    if (doGPUGrid) {
      myGrid = (ObitUVGrid*)outImage->myGrid;
      ifacet = 0; nfacet=1;
      /* Maybe not if (doBeam) {nfacet = 2; ifacet = 0;}*/
      if (!myGrid->gridInfo)
	myGrid->gridInfo = ObitGPUGridCreate("GPUGrid", nfacet, num_GPU, cuda_device, 
					     (Obit*)outImage, inUV, TRUE);
      iGPU = myGrid->gridInfo->cudaInfo->FacetGPU[ifacet]; /* which GPU to do this in? */
      ObitGPUGridSetGPUStruct(myGrid->gridInfo, (Obit*)myGrid, chDone, ifacet, nfacet, 
			      iGPU, num_GPU, cuda_device[iGPU], nplane, inUV, (Obit*)outImage, 
			      doBeam, err);
      if (err->error) Obit_traceback_msg (err, routine, outImage->name);
      myGrid->gridInfo->cudaInfo->gpu_memory[iGPU] = (size_t)gpumem[iGPU];
      /*???myGrid->gridInfo->GPU_device_no[iGPU] = cuda_device[iGPU];*/
      /*??myGrid->gridInfo->cudaInfo->GPU_device_no[iGPU] = cuda_device[iGPU];*/
    }

  /* Reset UV buffer size to at least 2 MByte per thread - not if GPU gridding */
  oldNPIO = 1024;  /* Save old value */
  ObitInfoListGetTest (inUV->info, "nVisPIO", &type, dim,  (gpointer)&oldNPIO);
  dim[0] = dim[1] = dim[2] = 1;
  if (!doGPUGrid) {
    /* How many threads? */
    myGrid = (ObitUVGrid*)outImage->myGrid;
    myGrid->nThreads = MAX (1, ObitThreadNumProc(outImage->thread));
    NPIO = myGrid->nThreads * 
      (olong) (0.5 + MAX (2.0, 2.0/(inUV->myDesc->lrec*4.0/(1024.0*1024.0))));
    ObitInfoListAlwaysPut (inUV->info, "nVisPIO",  OBIT_long, dim,  (gpointer)&NPIO);
   /* end not doGPUGrid */
  } else {  /* GPUGrid */
    NPIO = GPU_NVISPIO;
    ObitInfoListAlwaysPut (inUV->info, "nVisPIO",  OBIT_long, dim,  (gpointer)&NPIO);
  } /* end doGPUGrid */
 
  /* close to force reset */
  ObitUVClose (inUV, err);
  /* Calibrating or selecting? */
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;
  ObitUVOpen (inUV, access, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  
  
  /* Grid Image */
  gridClass->ObitUVGridReadUV ((ObitUVGrid*)outImage->myGrid, inUV, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  
  /* Gridding correction - write */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (outImage->info, "Channel",  OBIT_long, dim,  &ichannel);
  gridClass->ObitUVGridFFT2Im((ObitUVGrid*)outImage->myGrid, (Obit*)outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  
  /* tell Max/Min */
  imMax = outImage->myDesc->maxval;
  imMin = outImage->myDesc->minval;
  Obit_log_error(err, OBIT_InfoErr, 
		 "Image max %g, min %g for %s", imMax, imMin,  outImage->name);
  /* ObitErrLog(err);  Progress Report */
    
  /* Copy Gaussian beam fit */
  if (theBeam->myDesc->beamMaj>0.0) {
    outImage->myDesc->beamMaj = theBeam->myDesc->beamMaj;
    outImage->myDesc->beamMin = theBeam->myDesc->beamMin;
    outImage->myDesc->beamPA  = theBeam->myDesc->beamPA;
  }
  
  /* Save last beam normalization in Beam infoList as "SUMWTS" */
  if (doBeam) {
    beamGrid = (ObitUVGrid*)(theBeam->myGrid);
    sumwts = beamGrid->BeamNorm;
    dim[0] = 1;
    dim[1] = 0;
    ObitInfoListPut(theBeam->info, "SUMWTS", OBIT_float, dim, (gpointer)&sumwts, err);
    if (err->error) Obit_traceback_msg (err, routine, theBeam->name);
  }
  
  /* Make sure Stokes correct */
  outImage->myDesc->crval[outImage->myDesc->jlocs] = 
    inUV->myDesc->crval[inUV->myDesc->jlocs];
  
  /* Close Image */
  ObitImageClose (outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  
  } /* end loop over channels */

  /* Free image buffer */
  outImage->image = ObitFArrayUnref(outImage->image);
  
  /* Free gridding members */
  outImage->myGrid = (Obit*)ObitUVGridUnref(outImage->myGrid);
  theBeam = imgClass->ObitImageGetBeam(outImage, iorder, plane, err);
  theBeam->myGrid  =(Obit*) ObitUVGridUnref(theBeam->myGrid);
  /* were the two the same? */
  if (theBeam->myGrid==NULL) outImage->myGrid = NULL;
  
  /* Reset NPIO */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (inUV->info, "nVisPIO",  OBIT_long, dim,  (gpointer)&oldNPIO);

  /* GPU arrays */
  if (cuda_device) g_free(cuda_device);
  if (gpumem)      g_free(gpumem);
}  /* end ObitImageUtilMakeImage */

/**
 * Parallel Grids, FFTs and makes corrections for the gridding convolution
 * for a set of images with a single read of the uv data.
 * Uses (creating if necessary) the myGrid member of out.
 * If calibration is being applied, a UV data opject will be open for each
 * image and beam; this may blow the limit on open file descriptors.
 * Can make multichannel continuum or spectral line images
 * \param inUV     Input uv data. Should be in form of stokes to be imaged
 *                 will all calibration and selection applied and any 
 *                 weighting applied.  If calibration is specified then
 *                 infoList element "ParGainUse" (int), if present gives image
 *                 specifies calibration table numbers.  
 *                 One value is used for all or oner per image
 * \param nPar     Number of parallel images
 * \param outImage Array of Images to be written.  Must be previously instantiated.
 *                 Beam normalization factor is written to output Beam
 *                 infoList as "SUMWTS"
 *                 Supports ObitImage, ObitImageMF or ObitImageWB with all the same order.
 * \param doBeam   if TRUE also make beam.  Will make the myBeam member of 
 *                 outImage.
 *                 If FALSE, and myGrid->BeamNorm 0.0 then reads SUMWTS value 
 *                 from beam infolist
 *                 If an image does not have a beam one is always made.
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
 * \li "do3D"   OBIT_bool (1,1,1) 3D image, else 2D? [def TRUE]
 * \li "doGPUGrid" OBIT_bool (1,1,1) use GPU for Gridding [def FALSE]
 * \li "chDone" OBIT_bool (nPar,1,1) Array of flags indicating channels are done [all FALSE]
 * \param err      Error stack, returns if not empty.
 */
void ObitImageUtilMakeImagePar (ObitUV *inUV, olong nPar, ObitImage **outImage, 
				gboolean doBeam, gboolean doWeight,
				ObitErr *err)
{
  ObitImage *theBeam=NULL;
  ObitImageDesc *IODesc;
  ofloat sumwts, imMax, imMin, BeamTaper=0.0, Beam[3]={0.0,0.0,0.0};
  gchar outName[120];
  ollong *gpumem=NULL;
  olong ldevice[20], *cuda_device=NULL, iGPU=-1, num_GPU=1; /* default GPU number*/
  olong i, j, ip, pln, ifacet, nfacet, nplane=1, nImage, nGain=0, *gainUse=NULL, gain;
  olong plane[5], NPIO, oldNPIO, *arrayNx=NULL, *arrayNy=NULL;
  olong doCalib = 0, norder, iorder, nadd=0, bmInc=1;
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gboolean forceBeam, doCalSelect=FALSE, initGPU=TRUE, *chDone=NULL;
  ObitIOAccess access;
  ObitUVGrid **grids=NULL;
  ObitUVGrid *myGrid=NULL;
  ObitImage **imArray=NULL;
  ObitUVGridClassInfo *gridClass;
  ObitImageClassInfo *imgClass;
  gboolean doGPUGrid=FALSE, doneBeam=FALSE;
  union ObitInfoListEquiv InfoReal; 
  gchar *routine = "ObitImageUtilMakeImagePar";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));

  /* Apply uniform weighting? */
  if (doWeight) ObitUVWeightData (inUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  
  imgClass  = (ObitImageClassInfo*)outImage[0]->ClassInfo;          /* Image class */
  
   /* Use GPU? */
  InfoReal.itg = (olong)doGPUGrid; type = OBIT_bool;
  ObitInfoListGetTest(inUV->info, "doGPUGrid", &type, (gint32*)dim, &InfoReal);
  doGPUGrid = InfoReal.itg;
  if (doGPUGrid) {
    /* How many GPUs? */
    if (ObitInfoListGetTest(inUV->info, "GPU_no", &type, (gint32*)dim, ldevice)) {
      /* Count, ignore -1s  */
      num_GPU = 0;
      for (i=0; i<dim[0]; i++) if (ldevice[i]>=0) num_GPU++;
      if (num_GPU<=0) {num_GPU=1; ldevice[0]=0;}  /* defaults to cude device 0 */
      cuda_device = g_malloc0((num_GPU+10)*sizeof(olong));
      gpumem =      g_malloc0((num_GPU+10)*sizeof(olong));
      for (i=0; i<num_GPU; i++) {
	if ((ldevice[i]<0) || (ldevice[i]>50)) ldevice[i] = i;  /* Sanity check */
	cuda_device[i] = ldevice[i];
	if (ObitGPUGridCheckGPU((olong)cuda_device[i])<0) { /* Check device number */
	  Obit_log_error(err, OBIT_Error, "Invalid GPU device %d %d",cuda_device[i],i);
	  return;
	} /* end check device */
	gpumem[i] =      ObitGPUGridSetGPU (cuda_device[i]); /* initialize */
      }
    } else {  /* default GPU = 0 */
      num_GPU = 1;
      cuda_device = g_malloc0((num_GPU+10)*sizeof(olong));
      gpumem      = g_malloc0((num_GPU+10)*sizeof(olong));
      cuda_device[0] = 0;
      if (ObitGPUGridCheckGPU((olong)cuda_device[0])<0) { /* Check device number */
	Obit_log_error(err, OBIT_Error, "%s: Invalid GPU device %d",routine,cuda_device[0]);
	return;
      } /* end check device */
      gpumem[0]      = ObitGPUGridSetGPU (cuda_device[0]); /* initialize */
    }
    if (err->prtLv>=2)
      Obit_log_error(err, OBIT_InfoErr, "Doing GPU Gridding  with %d GPUs",num_GPU);
  } /* end using GPU */

  /* Need new gridding member? */
  /* How many to do */
  norder = imgClass->ObitImageGetBeamOrder(outImage[0]); /* No. orders of beam */
  nImage = nPar;
  if (doBeam) {
    nImage += (norder+1)*nImage;  /* Doing Beams as well? */
    bmInc = 2;  /* increment in grids array */
  } else {
    /* Add number of images without a beam  */
    for (j=0; j<nPar; j++) {
      sumwts = 0.0;
      theBeam = (ObitImage*)outImage[j]->myBeam;
      if ((theBeam==NULL) ||
	  !ObitInfoListGetTest(theBeam->info, "SUMWTS", &type, dim, &sumwts)) {
	nImage += (norder+1);
	nadd++;  /* How many beams being added */
      }
    }
  }
  
  /* Calibration */
  doCalib = 0;
  ObitInfoListGetTest(inUV->info, "doCalib", &type, dim, &doCalib);
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  dim[0] = 0;
  ObitInfoListGetP(inUV->info, "ParGainUse", &type, dim, (gpointer*)&gainUse);
  nGain = dim[0];
  /* Ch done array? */
  ObitInfoListGetP(inUV->info, "chDone", &type, dim, (gpointer*)&chDone);
  
  /* Create work arrays */
  grids   = g_malloc0(nImage*sizeof(ObitUVGrid*));
  imArray = g_malloc0(nImage*sizeof(ObitImage*));
  arrayNx = g_malloc0(nImage*sizeof(olong));
  arrayNy = g_malloc0(nImage*sizeof(olong));
  
  /* Ensure data fully instantiated and OK */
  ObitUVFullInstantiate (inUV, TRUE, err);
  if (err->error) goto cleanup;
  
  /* Loop over images creating things */
  initGPU = TRUE;  /* Need to set up for using GPU */
  ip = 0; ifacet = 0;
  for (j=0; j<nPar; j++) {
    
    /*  Open image ReadOnly to get proper descriptor, pixel array */
    dim[0] = 7;
    for (i=0; i<IM_MAXDIM; i++) {blc[i] = 1; trc[i] = 0;}
    ObitInfoListPut (outImage[j]->info, "BLC", OBIT_long, dim, blc, err); 
    ObitInfoListPut (outImage[j]->info, "TRC", OBIT_long, dim, trc, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage[j]->name);
    outImage[j]->image = ObitImageUnref(outImage[j]->image);  /* No buffer needed here */
    if ((ObitImageOpen (outImage[j], OBIT_IO_ReadOnly, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, 
		     "ERROR opening image %s", outImage[j]->name);
      goto cleanup;
    }
    ObitImageClose (outImage[j], err);
    
    /* Need new gridding member? */
    sprintf (outName, "UVGrid for %s ",outImage[j]->name);
    IODesc = (ObitImageDesc*)outImage[j]->myIO->myDesc;
    if (outImage[j]->myGrid == NULL) {
      if (ObitImageWBIsA(outImage[j])) {
	outImage[j]->myGrid = (Obit*)newObitUVGridWB(outName);
	nplane = 1;  /* Note GPU gridding doesn't work for this */
      } else if (ObitImageMFIsA(outImage[j])) {
	outImage[j]->myGrid = (Obit*)newObitUVGridMF(outName);
	nplane = ((ObitImageMF*)outImage[j])->nSpec;
     /* Spectroscopic Image cube? Use ImageMF gridder */
      } else if (IODesc->inaxes[IODesc->jlocf]>1) {
	outImage[j]->myGrid = (Obit*)newObitUVGridMF(outName);
	nplane = IODesc->inaxes[IODesc->jlocf]; 
      }  else {
	outImage[j]->myGrid = (Obit*)newObitUVGrid(outName);
	nplane = 1;
      }
    } /* end new grids */
    
    /* Get Beam member from outImage[j] */
    theBeam  = (ObitImage*)outImage[j]->myBeam;
    /* Need to force making a beam? */
    forceBeam = (theBeam==NULL) ||
      !ObitInfoListGetTest(theBeam->info, "SUMWTS", &type, dim, (gpointer)&sumwts);
    /* Arrays -  Need set for the Beam? */
    if (doBeam || forceBeam) {
      /* Ensure beam image OK */
      theBeam->image = ObitImageUnref(theBeam->image);  /* No buffer needed here */
      ObitImageOpen (theBeam, OBIT_IO_WriteOnly, err);
      ObitImageClose (theBeam, err);
      if (err->error) goto cleanup;
      
      /* Loop over beams */
      for (iorder=0; iorder<=norder; iorder++) {
	/*array[ip] = theBeam->image;*/
	imArray[ip] = ObitImageRef(theBeam);
	arrayNx[ip] = theBeam->myDesc->inaxes[0];
	arrayNy[ip] = theBeam->myDesc->inaxes[1];
	sprintf (outName, "UVGrid for %s Beam %d",outImage[j]->name,iorder);
	IODesc = (ObitImageDesc*)outImage[j]->myIO->myDesc;
	/* WB or MF Image? */
	if (ObitImageWBIsA(outImage[j])) {
	  grids[ip] = (ObitUVGrid*)newObitUVGridWB(outName);
	} else if (ObitImageMFIsA(outImage[j])) {
	  grids[ip] = (ObitUVGrid*)newObitUVGridMF(outName);
	  /* Spectroscopic Image cube? - use MF Gridding */
	} else if (IODesc->inaxes[IODesc->jlocf]>1) {
	    grids[ip] = (ObitUVGrid*)newObitUVGridMF(outName);
	} else {
	  grids[ip] = newObitUVGrid(outName);
	}
	/* Save which beam on gridder */
	pln = 1+iorder;
	dim[0] = dim[1] = dim[2] = 1;
	ObitInfoListAlwaysPut (grids[ip]->info, "Channel",  OBIT_long, dim,  &pln);
	/* Any calibration setup */
	if ((doCalib>0) && (ip==0) ) {
	  gain = 0;
	  if (nGain==1)    gain = gainUse[0];
	  if (nGain>=nPar) gain = gainUse[j];
	  dim[0] = 1; dim[1] = 1; dim[2] = 1;
	  ObitInfoListAlwaysPut(inUV->info, "gainUse", OBIT_long, dim, &gain);
	  ObitInfoListAlwaysPut(inUV->info, "doCalib", OBIT_long, dim, &doCalib);
	}
	/* Gridding setup */
	gridClass = (ObitUVGridClassInfo*)grids[ip]->ClassInfo; /* Gridder class */
	gridClass->ObitUVGridSetup (grids[ip], inUV, (Obit*)theBeam, TRUE, err);
	/* WB Image? Save beam order */
	if (ObitImageWBIsA(outImage[j])) ((ObitUVGridWB*)grids[ip])->order = iorder;
	if ((doCalib<=0) && (ip==0)) {  /* Need files open for calibration */
	  ObitUVClose(inUV, err);
	  if (err->error) goto cleanup;
	}
	/* This one done? -CHECK THIS 
	if (chDone) ((ObitUVGrid*)grids[ip])->isDone = chDone[j];*/
	ip++;
	if (err->error) goto cleanup;
      } /* end loop over beams (iorder) */
    } /* end if doBeam */
  
    /* Check that ip in bounds */
    Obit_return_if_fail((ip<nImage), err, "%s: Blew internal array %d >= %d", 
			routine, ip, nImage) ;
    /* Image grid */
    sprintf (outName, "UVGrid for %s image",outImage[j]->name);
    
    imArray[ip] = ObitImageRef(outImage[j]);
    pln = 1;
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (imArray[ip]->info, "Channel",  OBIT_long, dim, &pln);
    arrayNx[ip] = outImage[j]->myDesc->inaxes[0];
    arrayNy[ip] = outImage[j]->myDesc->inaxes[1];
    /*array[ip] = outImage[j]->image;*/
    grids[ip]   = ObitUVGridRef(outImage[j]->myGrid);
    /* Gridding setup */
    gridClass = (ObitUVGridClassInfo*)grids[ip]->ClassInfo; /* Gridder class */
    gridClass->ObitUVGridSetup (grids[ip], inUV, (Obit*)outImage[j], FALSE, err);
    if (doCalib<=0) {  /* Need files open for calibration - 
			  this may blow limit on open files */
      ObitUVClose(inUV, err);
      if (err->error) goto cleanup;
    }
    /* This one done? -CHECK THIS
    if (chDone) ((ObitUVGrid*)grids[ip])->isDone = chDone[j]; */
    ip++;
    if (err->error) goto cleanup;

    /* GPU Info on first grid */
    if (doGPUGrid) {
      iGPU = (iGPU+1)%num_GPU; /* which GPU to do this in? */
      myGrid = grids[0];  /* All CUDA info on grids[0]->gridInfo->cudaInfo */
      nfacet = nPar;
      if (doBeam) nfacet *= 2;
      else        nfacet += nadd;
      if (initGPU && (myGrid->gridInfo==NULL)) {  /* Init GPU on first facet */
	myGrid->gridInfo = ObitGPUGridCreate("GPUGrid", nfacet, num_GPU, cuda_device, 
					     (Obit*)outImage[0], inUV, TRUE);
	myGrid->doGPUGrid = doGPUGrid; 
	initGPU = FALSE;  /* Don't do again */
      }
      /* Also beam? Beam preceeds image */
      if (doBeam || forceBeam) {
	theBeam = (ObitImage*)outImage[j]->myBeam; 
	theBeam->myGrid = ObitUVGridRef(outImage[j]->myGrid);  /* Link grid info */
	if (!doneBeam) { /* Only one */
	  iGPU = myGrid->gridInfo->cudaInfo->FacetGPU[ifacet]; /* which GPU to do this in? */
	  ObitGPUGridSetGPUStruct (myGrid->gridInfo, (Obit*)grids[j*bmInc], chDone,
				   ifacet, nfacet, iGPU, num_GPU, cuda_device[iGPU], nplane, inUV, 
				   (Obit*)theBeam, TRUE, err);
	  //doneBeam = TRUE;
	}
	if (err->error) Obit_traceback_msg (err, routine, outImage[j]->name);
	theBeam->myGrid = ObitUVGridUnref(theBeam->myGrid);
	ifacet++;
      }
      /* Image */
      iGPU = myGrid->gridInfo->cudaInfo->FacetGPU[ifacet]; /* which GPU to do this in? */
      ObitGPUGridSetGPUStruct(myGrid->gridInfo, (Obit*)grids[j*bmInc],  chDone,
			      ifacet, nfacet, iGPU, num_GPU, cuda_device[iGPU], nplane, inUV, 
			      (Obit*)outImage[j], FALSE, err);
      ifacet++;
      if (err->error) Obit_traceback_msg (err, routine, outImage[j]->name);
      myGrid->gridInfo->cudaInfo->cuda_device[iGPU] = cuda_device[iGPU];
      myGrid->gridInfo->cudaInfo->gpu_memory[iGPU] = (size_t)gpumem[iGPU];
    } /* end doGPUGrid */
    
  } /* end loop over images */
  
  /* get any previous normalization  */
  for (j=0; j<nPar; j++) {
    sumwts = 0.0;
    myGrid = (ObitUVGrid*)outImage[j]->myGrid;
    if (doBeam) {
      theBeam = (ObitImage*)outImage[j]->myBeam;
      if ((theBeam!=NULL) &&
	  ObitInfoListGetTest(theBeam->info, "SUMWTS", &type, dim, &sumwts)) {
	/* Save beam normalization in Beam */
	myGrid->BeamNorm = sumwts;
      }
    } else { /* Reset max/min values */
      /* Open to update header with default max, min */
      theBeam = (ObitImage*)outImage[j]->myBeam;
      theBeam->extBuffer = TRUE;   /* No need for buffer */
      ObitImageOpen (theBeam, OBIT_IO_ReadWrite, err);
    if (err->error) goto cleanup;
      theBeam->myDesc->maxval = -1.0e20;
      theBeam->myDesc->minval =  1.0e20;
      ((ObitImageDesc*)theBeam->myIO->myDesc)->maxval = -1.0e20;
      ((ObitImageDesc*)theBeam->myIO->myDesc)->minval =  1.0e20;
      theBeam->myStatus = OBIT_Modified;  /* force update */
      ObitImageClose (theBeam, err);
      if (err->error) goto cleanup;
      theBeam->extBuffer = FALSE;   /* May need for buffer */
    }
    /* Open to update header with default max, min */
    outImage[j]->extBuffer = TRUE;   /* No need for buffer */
    ObitImageOpen (outImage[j], OBIT_IO_ReadWrite, err);
    if (err->error) goto cleanup;
    /* Reset max/min values */
    outImage[j]->myDesc->maxval = -1.0e20;
    outImage[j]->myDesc->minval =  1.0e20;
    ((ObitImageDesc*)outImage[j]->myIO->myDesc)->maxval = -1.0e20;
    ((ObitImageDesc*)outImage[j]->myIO->myDesc)->minval =  1.0e20;
    outImage[j]->myStatus = OBIT_Modified;  /* force update */
    ObitImageClose (outImage[j], err);
    if (err->error) goto cleanup;
    outImage[j]->extBuffer = FALSE;   /* May need for buffer */
  }  /* end loop over images */
  
  /* Applying calibration or selection? */
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;
  
  /* Reset UV buffer size to at least 2 MByte per thread */
  oldNPIO  = 1024;  /* Save old value */
  ObitInfoListGetTest (inUV->info, "nVisPIO", &type, dim,  (gpointer)&oldNPIO);
  dim[0] = dim[1] = dim[2] = 1;
  if (!doGPUGrid) {
    /* How many threads? */
    grids[0]->nThreads = MAX (1, ObitThreadNumProc(grids[0]->thread));
    NPIO = ObitImageUtilBufSize(inUV);
    ObitInfoListAlwaysPut (inUV->info, "nVisPIO",  OBIT_long, dim,  (gpointer)&NPIO);
  } else {  /* GPUGrid */
    NPIO = GPU_NVISPIO;
    ObitInfoListAlwaysPut (inUV->info, "nVisPIO",  OBIT_long, dim,  (gpointer)&NPIO);
    /* Close to make sure updated */
    ObitUVClose(inUV, err);
    if (err->error) goto cleanup;
 } /* end doGPUGrid */
  
  /*  Open uv data for first if needed */
  if ((inUV->myStatus!=OBIT_Active) && (inUV->myStatus!=OBIT_Modified)) {
    ObitUVOpen (inUV, access, err);
    if (err->error) goto cleanup;
  }
  
  /* Make images 
     Grid  */
  gridClass = (ObitUVGridClassInfo*)grids[0]->ClassInfo; /* Gridder class */
  if (nImage>1) 
    gridClass->ObitUVGridReadUVPar (nImage, grids, inUV, err);
  else /* Only one */
    gridClass->ObitUVGridReadUV (grids[0], inUV, err);
  if (err->error) goto cleanup;
  
  /* Close uv data files if needed */
  if (doCalib>0) {  /* Need files open for calibration */
    for (ip=0; ip<nImage; ip++) {
      ObitUVClose(inUV, err);
      if (err->error) goto cleanup;
    }
  }
  
  /* In case beam normalization lost... */
  for (j=0; j<nPar; j++) {
    theBeam = (ObitImage*)outImage[j]->myBeam;
    ObitInfoListGetTest(theBeam->info, "SUMWTS", &type, dim, (gpointer)&sumwts);
    if (grids[j]->BeamNorm<=1.0)  grids[j]->BeamNorm = sumwts;
  }
    /* FFT, Gridding correction - deal with beam normalization here 
     Write output images */
  gridClass->ObitUVGridFFT2ImPar (nImage, grids, (Obit**)imArray, err);
  if (err->error) goto cleanup;
  
  /* Get restoring beam */
  ObitInfoListGetTest(inUV->info, "Beam", &type, dim, Beam);

  /* Loop over images finishing */
  ip = 0;
  for (j=0; j<nPar; j++) {
    
    /* Get Beam member from outImage[j], test if already made */
    theBeam = (ObitImage*)outImage[j]->myBeam;
    forceBeam = (theBeam==NULL) ||
      !ObitInfoListGetTest(theBeam->info, "SUMWTS", &type, dim, (gpointer)&sumwts);
    /* In case beam normalization lost... */
    if (grids[j]->BeamNorm<=1.0)  grids[j]->BeamNorm = sumwts;
    
    /* Next a beam? */
    if (doBeam || forceBeam) {
      /* Loop over beams */
      for (iorder=0; iorder<=norder; iorder++) {
	/*if (((ObitUVGrid*)grids[j])->isDone) continue;   Already finished? */
	/* Get relevant beam */
	theBeam = imgClass->ObitImageGetBeam(outImage[j], iorder, plane, err);
	if (err->error) goto cleanup;
	
	if (iorder==0) {  /* Stuff for dirty beam */
	  /* Tell Sum of gridding weights */
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Sum of Weights %g for %s",
			 grids[j]->BeamNorm,outImage[j]->name);
	  
	  /* tell Max/Min */
	  imMax = theBeam->myDesc->maxval;
	  imMin = theBeam->myDesc->minval;
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Beam max %g, min %g for %s", imMax, imMin, theBeam->name);
	  
	  /* Fit beam */
	  ObitImageUtilFitBeam (theBeam, err);
	  if (err->error) goto cleanup;

	  /* Beam specified? */
	  ObitInfoListGetTest(outImage[j]->myDesc->info, "BeamTapr", &type, dim, &BeamTaper);
	  if ((Beam[0]>0.0) && (BeamTaper<=0.0)) {
	     if (j==0) Obit_log_error(err, OBIT_InfoErr, 
	     "Using Beam %f %f %f", Beam[0], Beam[1], Beam[2]);
	     theBeam->myDesc->beamMaj = Beam[0]/3600.0;
	     theBeam->myDesc->beamMin = Beam[1]/3600.0;
	     theBeam->myDesc->beamPA  = Beam[2];
	  } else if (Beam[0]<=0.0) {  /* Use first fitted beam if none given */
	    Beam[0] = theBeam->myDesc->beamMaj*3600.0;
	    Beam[1] = theBeam->myDesc->beamMin*3600.0;
	    Beam[2] = theBeam->myDesc->beamPA;
	    /* Save */
	    dim[0] = 3; dim[1]= dim[2] = dim[3] = dim[4] = 1;
	    ObitInfoListAlwaysPut(inUV->info, "Beam", OBIT_float, dim, Beam);
	  }
	  
	  /* Save last beam normalization in Beam infoList as "SUMWTS" */
	  myGrid = (ObitUVGrid*)outImage[j]->myGrid;
	  sumwts = myGrid->BeamNorm;
	  dim[0] = 1; dim[1] = 1;
	  ObitInfoListAlwaysPut(theBeam->info, "SUMWTS", OBIT_float, dim, (gpointer)&sumwts);
	} /* End dirty beam stuff */
	
	/* Free image buffer */
	theBeam->image = ObitImageUnref(theBeam->image);
	imArray[ip] = ObitImageUnref(imArray[ip]);  /* Free image */
	
	ip++;  /* update image array pointer */
      } /* end loop over beams */
    } /* end if doBeam */
    
    /* if (!((ObitUVGrid*)grids[ip])->isDone) { Already finished? */
      /* Now image - open to update header with restoring beam */
      ObitImageOpen (outImage[j], OBIT_IO_ReadWrite, err);
      if (err->error) goto cleanup;
      outImage[j]->myStatus = OBIT_Modified;  /* force update */
      
      /* tell Max/Min */
      imMax = imArray[ip]->myDesc->maxval;
      imMin = imArray[ip]->myDesc->minval;
      Obit_log_error(err, OBIT_InfoErr, 
		     "Image max %g, min %g for %s", imMax, imMin,  outImage[j]->name);
      ObitErrLog(err); 
      
      /* Copy Gaussian beam fit if necessary */
      theBeam = (ObitImage*)outImage[j]->myBeam;
      /*if ((outImage[j]->myDesc->beamMaj<=0.0) || (outImage[j]->myDesc->beamMin<=0.0)) {*/
      if (theBeam->myDesc->beamMaj>0.0) {
	outImage[j]->myDesc->beamMaj = theBeam->myDesc->beamMaj;
	outImage[j]->myDesc->beamMin = theBeam->myDesc->beamMin;
	outImage[j]->myDesc->beamPA  = theBeam->myDesc->beamPA;
      }
      
      /* Make sure Stokes correct if one of I,Q,U,V */
      if (inUV->myDesc->crval[inUV->myDesc->jlocs]>0.0)
	outImage[j]->myDesc->crval[outImage[j]->myDesc->jlocs] = 
	  inUV->myDesc->crval[inUV->myDesc->jlocs];
      
      ObitImageClose (outImage[j], err);
      if (err->error) goto cleanup;
      imArray[ip] = ObitImageUnref(imArray[ip]);  /* Free image */
      /* }  end not already finished */
    ip++;
  } /* end loop writing output */
  
    /* Cleanup */
 cleanup:
  
  /* Reset NPIO */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (inUV->info, "nVisPIO",  OBIT_long, dim,  (gpointer)&oldNPIO);

  /* Free gridding objects */
  if (grids) {
    for (j=0; j<nImage; j++) {
      grids[j] = ObitUVGridUnref(grids[j]);
    }
    g_free(grids); grids = NULL;
  }
  for (j=0; j<nPar; j++) {
    outImage[j]->myGrid = ObitUVGridUnref(outImage[j]->myGrid);
  }
  /* Free Image Buffers */
  for (j=0; j<nPar; j++) {
    outImage[j]->image = ObitFArrayUnref(outImage[j]->image);
    if (outImage[j]->myBeam)
      ((ObitImage*)outImage[j]->myBeam)->image = 
	ObitFArrayUnref(((ObitImage*)outImage[j]->myBeam)->image);
   }
 
  g_free(imArray); 
  g_free(arrayNx);
  g_free(arrayNy);
  /* GPU arrays */
  if (cuda_device) g_free(cuda_device);
  if (gpumem)      g_free(gpumem);
  
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
}  /* end ObitImageUtilMakeImagePar */

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
  olong iblc[IM_MAXDIM], itrc[IM_MAXDIM], oblc[IM_MAXDIM], otrc[IM_MAXDIM];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, j;
  odouble RAPnt, DecPnt;
  olong nTh, nrow, lorow, hirow, nrowPerThread, nThreads;
  InterpFuncArg **threadArgs;
  ofloat fblank = ObitMagicF();
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
    /* Index as well */
    ObitImageDescIndex(outImage->myDesc);

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
  /* Index to be sure */
  ObitImageDescIndex(outImage->myDesc);

  /* Read input plane */
  if ((ObitImageRead (inImage,NULL , err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR reading image %s", 
		   routine, inImage->name);
    return;
  }

  /* Convert pure zero to fblank */
  ObitFArrayInClip (inImage->image, -1.0e-25, 1.0e-25, fblank);

  /* Make interpolator */
  interp = newObitFInterpolateCreate ("Interpolator", inImage->image, 
				      inImage->myDesc, hwidth);

   /* Initialize Threading */
  nThreads = MakeInterpFuncArgs (inImage->thread, -1, 0, NULL,
				 inImage->myDesc, outImage->myDesc, 
				 interp, err, &threadArgs);

  /* Divide up work */
  nrow = outImage->myDesc->inaxes[1];
  nrowPerThread = nrow/nThreads;
  nTh = nThreads;
  if (nrow<64) {nrowPerThread = nrow; nTh = 1;}
  /* No fewer than 64 rows per thread */
  if ((nrowPerThread)<64) {
    nTh = (olong)(0.5 + nrow/64.0);
   if (nTh>0)  nrowPerThread = nrow/nTh;
  }
  if (nTh<=0) {nrowPerThread = nrow; nTh = 1;}
  lorow = 1;
  hirow = nrowPerThread;
  hirow = MIN (hirow, nrow);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hirow = nrow;  /* Make sure do all */
    threadArgs[i]->inData  = ObitFArrayRef(inImage->image);
    threadArgs[i]->outData = ObitFArrayRef(outImage->image);
    threadArgs[i]->first   = lorow;
    threadArgs[i]->last    = hirow;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which row */
    lorow += nrowPerThread;
    hirow += nrowPerThread;
    hirow = MIN (hirow, nrow);
  }

  /* Do operation */
  OK = ObitThreadIterator (inImage->thread, nTh, 
			   (ObitThreadFunc)ThreadImageInterp,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
  /* */

  /* Free local objects */
  interp = ObitFInterpolateUnref(interp);
  KillInterpFuncArgs(nThreads, threadArgs);

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
 * in inImage.
 * There is no interpolation between planes
 * \param inImage  Image to be interpolated.
 * \param outImage Image to be written.  Must be previously instantiated.
 * \param XPix     Image of x pixels in inImage for outImage
 * \param YPix     Image of y pixels in inImage for outImage
 * \param inPlane  desired plane in inImage, 1-rel pixel numbers on planes 3-7
 * \param outPlane desired plane in outImage
 * \param hwidth   interpolation halfwidth (1 or 2 usually OK, 4 max)
 * \param err      Error stack, returns if not empty.
 */
void 
ObitImageUtilInterpolateImageXY (ObitImage *inImage, ObitImage *outImage, 
				 ObitImage *XPix, ObitImage *YPix, 
			         olong *inPlane, olong *outPlane,
			         olong hwidth, ObitErr *err)
{
  ObitIOSize IOBy;
  ObitFInterpolate *interp=NULL;
  ObitImageDesc *tmpDesc=NULL;
  olong iblc[IM_MAXDIM], itrc[IM_MAXDIM], oblc[IM_MAXDIM], otrc[IM_MAXDIM];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, j, pln[5];
  odouble RAPnt, DecPnt;
  olong nTh, nrow, lorow, hirow, nrowPerThread, nThreads;
  InterpFuncArg **threadArgs;
  ofloat fblank = ObitMagicF();
  gboolean OK;
  gchar *today=NULL;
  gchar *routine = "ObitImageUtilInterpolateImageXY";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));
  g_assert (inPlane!=NULL);
  g_assert (outPlane!=NULL);

  /* Precomputed pixels */
  /* Open images */
  if ((ObitImageOpen (XPix, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, inImage->name);
    return;
  }
  if ((ObitImageOpen (YPix, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, outImage->name);
    return;
  }
  /* Read */
  pln[0] = pln[1] = pln[2] = pln[3] = pln[4] = 1;
  ObitImageGetPlane (XPix, NULL, pln, err);
  if (err->error) Obit_traceback_msg (err, routine, XPix->name);
  ObitImageGetPlane (YPix, NULL, pln, err);
  if (err->error) Obit_traceback_msg (err, routine, YPix->name);
 
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
    /* Index as well */
    ObitImageDescIndex(outImage->myDesc);

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

  /* Make interpolator */
  interp = newObitFInterpolateCreate ("Interpolator", inImage->image, 
				      inImage->myDesc, hwidth);

   /* Initialize Threading */
  nThreads = MakeInterpFuncArgs (inImage->thread, -1, 0, NULL,
				 inImage->myDesc, outImage->myDesc, 
				 interp, err, &threadArgs);

  /* Divide up work */
  nrow = outImage->myDesc->inaxes[1];
  nrowPerThread = nrow/nThreads;
  nTh = nThreads;
  if (nrow<64) {nrowPerThread = nrow; nTh = 1;}
  /* No fewer than 64 rows per thread */
  if ((nrowPerThread)<64) {
    nTh = (olong)(0.5 + nrow/64.0);
   if (nTh>0)  nrowPerThread = nrow/nTh;
  }
  if (nTh<=0) {nrowPerThread = nrow; nTh = 1;}
  lorow = 1;
  hirow = nrowPerThread;
  hirow = MIN (hirow, nrow);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hirow = nrow;  /* Make sure do all */
    threadArgs[i]->inData  = ObitFArrayRef(inImage->image);
    threadArgs[i]->outData = ObitFArrayRef(outImage->image);
    threadArgs[i]->XPixData = ObitFArrayRef(XPix->image);
    threadArgs[i]->YPixData = ObitFArrayRef(YPix->image);
    threadArgs[i]->first   = lorow;
    threadArgs[i]->last    = hirow;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which row */
    lorow += nrowPerThread;
    hirow += nrowPerThread;
    hirow = MIN (hirow, nrow);
  }

  /* Do operation */
  OK = ObitThreadIterator (inImage->thread, nTh, 
			   (ObitThreadFunc)ThreadImageInterp,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
  /* */

  /* Free local objects */
  interp = ObitFInterpolateUnref(interp);
  KillInterpFuncArgs(nThreads, threadArgs);

  /* Write output */
  if ((ObitImageWrite (outImage, NULL, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR writing image %s", 
		   routine, outImage->name);
    return;
  }

  /* Close */
  if ((ObitImageClose (XPix, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, XPix->name);
    return;
  }
  if ((ObitImageClose (YPix, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, YPix->name);
    return;
  }
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
} /* end  ObitImageUtilInterpolateImageXY */

/**
 * Get input pixels in InImage for outImage.
 * \param inImage  Image to be interpolated.
 * \param outImage Image to be written.
 * \param XPix     Image of x pixels in inImage for outImage
 * \param YPix     Image of y pixels in inImage for outImage
 * \param err      Error stack, returns if not empty.
 */
void 
ObitImageUtilGetXYPixels (ObitImage *inImage, ObitImage *outImage, 
			  ObitImage *XPix, ObitImage *YPix, 
			  ObitErr *err)
{
  ObitImageDesc *tmpDesc=NULL;
  olong i, j;
  odouble RAPnt, DecPnt;
  olong nTh, nrow, lorow, hirow, nrowPerThread, nThreads;
  InterpFuncArg **threadArgs;
  gboolean OK;
  gchar *today=NULL;
  gchar *routine = "ObitImageUtilGetXYPixels";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));
  /* Are outImage and  ?Pix are compatable */
  Obit_return_if_fail((ObitFArrayIsCompatable(outImage->image, XPix->image)), err, 
		      "%s: output and XPix images incompatable", routine);
  Obit_return_if_fail((ObitFArrayIsCompatable(outImage->image, YPix->image)), err, 
		      "%s: output and YPix images incompatable", routine);

 
  /* Open images */
  if ((ObitImageOpen (XPix, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, inImage->name);
    return;
  }
  if ((ObitImageOpen (YPix, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, outImage->name);
    return;
  }
  /* Adjust out put descriptor  - copy from input */
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

  /* Index as well */
  ObitImageDescIndex(outImage->myDesc);

  /* Initialize Threading */
  nThreads = MakeInterpFuncArgs (inImage->thread, -1, 0, NULL,
				 inImage->myDesc, outImage->myDesc, 
				 NULL, err, &threadArgs);

  /* Divide up work */
  nrow = outImage->myDesc->inaxes[1];
  nrowPerThread = nrow/nThreads;
  nTh = nThreads;
  if (nrow<64) {nrowPerThread = nrow; nTh = 1;}
  /* No fewer than 64 rows per thread */
  if ((nrowPerThread)<64) {
    nTh = (olong)(0.5 + nrow/64.0);
   if (nTh>0)  nrowPerThread = nrow/nTh;
  }
  if (nTh<=0) {nrowPerThread = nrow; nTh = 1;}
  lorow = 1;
  hirow = nrowPerThread;
  hirow = MIN (hirow, nrow);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hirow = nrow;  /* Make sure do all */
    threadArgs[i]->inData  = ObitFArrayRef(inImage->image);
    threadArgs[i]->outData = ObitFArrayRef(outImage->image);
    threadArgs[i]->XPixData = ObitFArrayRef(XPix->image);
    threadArgs[i]->YPixData = ObitFArrayRef(YPix->image);
    threadArgs[i]->first   = lorow;
    threadArgs[i]->last    = hirow;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which row */
    lorow += nrowPerThread;
    hirow += nrowPerThread;
    hirow = MIN (hirow, nrow);
  }

  /* Do operation */
  OK = ObitThreadIterator (inImage->thread, nTh, 
			   (ObitThreadFunc)ThreadGetXYPixels,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);

  /* Free local objects */
  KillInterpFuncArgs(nThreads, threadArgs);

  /* Write outputs */
  if ((ObitImageWrite (XPix, NULL, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR writing image %s", 
		   routine, XPix->name);
    return;
  }
  if ((ObitImageWrite (YPix, NULL, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR writing image %s", 
		   routine, YPix->name);
    return;
  }

  /* Close */
  if ((ObitImageClose (XPix, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, XPix->name);
    return;
  }
  if ((ObitImageClose (YPix, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, YPix->name);
    return;
  }
  /* Free image buffer if not memory resident */
  if (XPix->mySel->FileType!=OBIT_IO_MEM) 
    XPix->image = ObitFArrayUnref(XPix->image);
  if (YPix->mySel->FileType!=OBIT_IO_MEM) 
    YPix->image = ObitFArrayUnref(YPix->image);

} /* end  ObitImageUtilGetXYPixels */

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
  olong iblc[IM_MAXDIM], itrc[IM_MAXDIM], oblc[IM_MAXDIM], otrc[IM_MAXDIM];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, j;
  gchar *today=NULL;
  olong nTh, nrow, lorow, hirow, nrowPerThread, nThreads;
  InterpFuncArg **threadArgs;
  gboolean OK;
  gchar *routine = "ObitImageUtilInterpolateImageZern";

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
    /* Index as well */
    ObitImageDescIndex(outImage->myDesc);

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
  
  /* Initialize Threading */
  nThreads = MakeInterpFuncArgs (inImage->thread, -1, nZern, ZCoef,
				 inImage->myDesc, outImage->myDesc, 
				 interp, err, &threadArgs);

  /* Divide up work */
  nrow = outImage->myDesc->inaxes[1];
  nrowPerThread = nrow/nThreads;
  nTh = nThreads;
  if (nrow<64) {nrowPerThread = nrow; nTh = 1;}
  /* No fewer than 64 rows per thread */
  if ((nrowPerThread)<64) {
    nTh = (olong)(0.5 + nrow/64.0);
    if (nTh>0) nrowPerThread = nrow/nTh;
  }
  if (nTh<=0) {nrowPerThread = nrow; nTh = 1;}
  lorow = 1;
  hirow = nrowPerThread;
  hirow = MIN (hirow, nrow);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hirow = nrow;  /* Make sure do all */
    threadArgs[i]->inData  = ObitFArrayRef(inImage->image);
    threadArgs[i]->outData = ObitFArrayRef(outImage->image);
    threadArgs[i]->first   = lorow;
    threadArgs[i]->last    = hirow;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which row */
    lorow += nrowPerThread;
    hirow += nrowPerThread;
    hirow = MIN (hirow, nrow);
  }

  /* Do operation */
  OK = ObitThreadIterator (inImage->thread, nTh, 
			   (ObitThreadFunc)ThreadImageInterp,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
  /* */

  /* Free local objects */
  interp = ObitFInterpolateUnref(interp);
  KillInterpFuncArgs(nThreads, threadArgs);

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
 * Also calculates a weight based on a circle defined by radiusfrom the center; 
 * this is 1.0 in the center and tapers with distance^2 to 0.0 outside.
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
				olong radius, olong *inPlane, olong *outPlane,
				olong hwidth, ObitErr *err)
{
  ObitIOSize IOBy;
  ObitFInterpolate *interp=NULL;
  olong blc[IM_MAXDIM], trc[IM_MAXDIM];
  gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong nTh, nrow, lorow, hirow, nrowPerThread, nThreads;
  InterpFuncArg **threadArgs;
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
 
  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;

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

  /* Check outImage and outWeight */
  if (!ObitFArrayIsCompatable(outImage->image, outWeight->image)) {
      Obit_log_error(err, OBIT_Error, "%s: Incompatable sizes for %s %s", 
		     routine, outImage->name, outWeight->name);
      return;
  }

  /* Make interpolator */
  interp = newObitFInterpolateCreate ("Interpolator", inImage->image, 
				      inImage->myDesc, hwidth);

  /* Initialize Threading */
  nThreads = MakeInterpFuncArgs (inImage->thread, radius, 0, NULL,
				 inImage->myDesc, outImage->myDesc, 
				 interp, err, &threadArgs);

  /* Divide up work */
  nrow = outImage->myDesc->inaxes[1];
  nrowPerThread = nrow/nThreads;
  nTh = nThreads;
  if (nrow<64) {nrowPerThread = nrow; nTh = 1;}
  /* No fewer than 64 rows per thread */
  if ((nrowPerThread)<64) {
    nTh = (olong)(0.5 + nrow/64.0);
    if (nTh>0) nrowPerThread = nrow/nTh;
  }
  if (nTh<=0) {nrowPerThread = nrow; nTh = 1;}
  lorow = 1;
  hirow = nrowPerThread;
  hirow = MIN (hirow, nrow);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hirow = nrow;  /* Make sure do all */
    threadArgs[i]->inData  = ObitFArrayRef(inImage->image);
    threadArgs[i]->outData = ObitFArrayRef(outImage->image);
    threadArgs[i]->wtData  = ObitFArrayRef(outWeight->image);
    threadArgs[i]->first   = lorow;
    threadArgs[i]->last    = hirow;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which row */
    lorow += nrowPerThread;
    hirow += nrowPerThread;
    hirow = MIN (hirow, nrow);
  }

  /* Do operation */
  OK = ObitThreadIterator (inImage->thread, nTh, 
			   (ObitThreadFunc)ThreadImageInterp,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
  /* */

  /* Free local objects */
  interp = ObitFInterpolateUnref(interp);
  KillInterpFuncArgs(nThreads, threadArgs);

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
} /* end  ObitImageUtilInterpolateWeight */

/**
 * Fill the pixels in outImage by corresponding values in inImage.
 * Like ObitImageUtilInterpolateWeight but without the interpolation.
 * Also calculates a weight based on a circle defined by radius from the center; 
 * this is 1.0 in the center and tapers with distance^2 to 0.0 outside.
 * If memOnly then the input image plane is assumed in inImage and only memory
 * resident parts of outImage and outWeight are modified.
 * \param inImage   Image to be copied
 * \param outImage  Image (*weight) to be written.  Must be previously instantiated.
 * \param outWeight Weight image to be written.  Must be previously instantiated and
 *                  have same geometry as outImage.
 * \param memOnly   if TRUE then work only in memory
 * \param radius    Radius in pixels of weighting circle
 * \param inPlane   Desired plane in inImage, 1-rel pixel numbers on planes 3-7; 
 *                  ignored if memOnly
 * \param outPlane  Desired plane in outImage; ignored if memOnly
 * \param err       Error stack, returns if not empty.
 * \return     TRUE if image needs to be interpolated, else FALSE.
 */
gboolean 
ObitImageUtilNoInterWeight (ObitImage *inImage, ObitImage *outImage, 
			    ObitImage *outWeight, gboolean memOnly,
			    olong radius, olong *inPlane, olong *outPlane,
			    ObitErr *err)
{
  gboolean doInter=TRUE;
  ObitIOSize IOBy;
  olong nx, ny, ix, iy, pos1[2], pos2[2], blc[IM_MAXDIM], trc[IM_MAXDIM];
  ollong indx;
  gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat inPixel[2], offPixel[2], xpos1[2], xpos2[2];
  ofloat *wtArr, dist2, rad2, irad2, xcen, ycen, wt;
  ofloat fblank = ObitMagicF();
  gboolean sameGrid;
  gchar *routine = "ObitImageUtilNoInterImage";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return doInter;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));

  /* Is interpolation needed? */
  if (!ObitImageDescAligned(inImage->myDesc, outImage->myDesc, err)) 
    return doInter;;
  if (err->error) Obit_traceback_val (err, routine, inImage->name, doInter);
  /* further sanity check, alignment pixel must be must be integer */
  inPixel[0] = outImage->image->naxis[0]/2; /* Use center of output */
  inPixel[1] = outImage->image->naxis[1]/2;
  ObitImageDescCvtPixel (outImage->myDesc, inImage->myDesc, inPixel, offPixel, err);
  if (err->error) Obit_traceback_val (err, routine, inImage->name, doInter);
  if (offPixel[0]>0.0) ix = (olong)(offPixel[0]+0.5);
  else                 ix = (olong)(offPixel[0]-0.5);
  sameGrid = fabs(offPixel[0]-ix)<0.001;
  if (offPixel[1]>0.0) iy = (olong)(offPixel[1]+0.5);
  else                 iy = (olong)(offPixel[1]-0.5);
  sameGrid = sameGrid && (fabs(offPixel[1]-iy)<0.001);
  if (!sameGrid) return doInter;
  
  /* Check outImage and outWeight */
  if (!ObitFArrayIsCompatable(outImage->image, outWeight->image)) {
      Obit_log_error(err, OBIT_Error, "%s: Incompatable sizes for %s %s", 
		     routine, outImage->name, outWeight->name);
      return doInter;
  }

  g_assert (ObitImageIsA(outWeight));
  g_assert (inPlane!=NULL);
  g_assert (outPlane!=NULL);
 
  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;

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
      return doInter;
    }
    if ((ObitImageOpen (outImage, OBIT_IO_ReadWrite, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		     routine, outImage->name);
      return doInter;
    }
    
    if ((ObitImageOpen (outWeight, OBIT_IO_ReadWrite, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		     routine, outWeight->name);
      return doInter;
    }
    /* Read input plane */
    if ((ObitImageRead (inImage,NULL , err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR reading image %s", 
		     routine, inImage->name);
      return doInter;
    }
  } /* end of not memory only */

  /* Do it -  init arrays */
  ObitFArrayFill(outImage->image,  0.0);
  ObitFArrayFill(outWeight->image, 0.0);

  /* Paste input into output array - need pixel alignment */
  pos1[0] = inImage->image->naxis[0]/2; /* Use center of input */
  pos1[1] = inImage->image->naxis[1]/2;
  xpos1[0] = pos1[0];  xpos1[1] = pos1[1];
  /* Find corresponding pixel in input */
  ObitImageDescCvtPixel (inImage->myDesc, outImage->myDesc, xpos1, xpos2, err);
  if (err->error) Obit_traceback_val (err, routine, inImage->name, doInter);
  pos2[0] = (olong)(xpos2[0]+0.5);  pos2[1] = (olong)(xpos2[1]+0.5);
 
  /* Copy image pixels */
  ObitFArrayShiftAdd (outImage->image, pos2, inImage->image, 
		      pos1, 1.0, outImage->image);
  /* Blank zeroes where there is no data*/
  ObitFArrayInClip(outImage->image, -1.0e-10, 1.0e-10, fblank); 

  /* Make Weight */
  /* Working version of radius */
  rad2 = radius * radius;
  irad2 = 1.0 / rad2;
  nx   = outWeight->myDesc->inaxes[0]; ny = outWeight->myDesc->inaxes[1];
  xcen = nx/2;  ycen = ny/2;

  /* Loop over image calculating weights */
  pos1[0] = pos1[1] = 0;
  wtArr = ObitFArrayIndex (outWeight->image, pos1);
  for (iy = 0; iy<ny; iy++) { /* loop in y */
    for (ix = 0; ix<nx; ix++) {/* loop in x */
      /* array index in out for this pixel */
      indx =  iy*nx + ix;
      /* weight based on distance from center */
      dist2 = (ix-xcen)*(ix-xcen) + (iy-ycen)*(iy-ycen);
      if (dist2 <= rad2) {
	wt = 1.0 - dist2 * irad2;
	wt = MAX (0.001, wt);
      } else wt = fblank;
      wtArr[indx] = wt;
    } /* end x array loop */
  } /* end y array loop */
 
  /* Multiply image by weight */
  ObitFArrayMul(outImage->image, outWeight->image, outImage->image);
  doInter = FALSE;  /* Not needed now */

  if (!memOnly) {
    /* Write image (*weight) output */
    if ((ObitImageWrite (outImage, NULL, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR writing image %s", 
		     routine, outImage->name);
      return doInter;
    }
    
    /* Close */
    /* Write weight */
    if ((ObitImageWrite (outWeight, NULL, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR writing image %s", 
		     routine, outWeight->name);
      return doInter;
    }
    
    if ((ObitImageClose (inImage, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		     routine, inImage->name);
      return doInter;
    }
    /* Free image buffer  if not memory resident */
    if (inImage->mySel->FileType!=OBIT_IO_MEM) 
      inImage->image = ObitFArrayUnref(inImage->image);
    if ((ObitImageClose (outImage, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		     routine, outImage->name);
      return doInter;
    }
    /* Free image buffer if not memory resident */
    if (outImage->mySel->FileType!=OBIT_IO_MEM) 
      outImage->image = ObitFArrayUnref(outImage->image);
    if ((ObitImageClose (outWeight, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		     routine, outImage->name);
      return doInter;
    }
    /* Free image buffer if not memory resident */
    if (outWeight->mySel->FileType!=OBIT_IO_MEM) 
      outWeight->image = ObitFArrayUnref(outWeight->image);
    
  } /* end of not memory only */

  return doInter;
} /* end ObitImageUtilNoInterWeight */

/**
 * Interpolate 3rd axis in inImage to inPlane and write to outImage outPlane
 * Input planes<1 or > nplanes will be blank filled.
 * \param inImage   Image to be interpolated. Honors BLC, TRC set
 * \param inPlane   fractional plane to interpolate (1-rel)
 * \param outImage  Image to be written.  Must be previously instantiated.
 * \param outPlane  plane to write in outImage (1-rel)
 * \param hwidth    Interpolation halfwidth (1 or 2 usually OK, 4 max)
 * \param err       Error stack, returns if not empty.
 */
void 
ObitImageUtilInterp3 (ObitImage *inImage, ofloat inPlane,
		      ObitImage *outImage, olong outPlane, 
		      olong hwidth, ObitErr *err)
{
  olong plane[5]={1,1,1,1,1}, i, nplane, first, npln;
  ofloat intWt[10];
  gchar *routine = "ObitImageUtilInterp3";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));

  /* Is inPlane in bounds? */
  npln = inImage->myDesc->inaxes[2];
  if ((inPlane<1.0) || (inPlane>npln)) {  /* Blank */
    if ((ObitImageOpen (outImage, OBIT_IO_ReadWrite, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		     routine, outImage->name);
      return;
    }
    /* blank */
    ObitFArrayFill(outImage->image, ObitMagicF());
    /* Write output */
    plane[0] = outPlane;
    if ((ObitImagePutPlane (outImage, NULL, plane, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR writing image %s", 
		     routine, outImage->name);
      return;
    }
    return;
  } /* end if output of bounds */
  
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
  /* Zero output FArray for accumulation */
  ObitFArrayFill(outImage->image, 0.0);
    
  /* Interpolation weights */
  SetConvKernal3 (inPlane, npln, hwidth, &first, intWt);

  nplane = 1 + hwidth*2;           /* How many planes to include */
  for (i=0; i<nplane; i++) {
    if (intWt[i]==0.0) continue;   /* don't bother if weight zero */
    /* Read input plane */
    plane[0] = first+i;
    if ((ObitImageGetPlane (inImage, NULL, plane, err) 
	 != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR reading image %s", 
		     routine, inImage->name);
      return;
    }
    /* Multiply input buffer by weights */
    ObitFArraySMul (inImage->image, intWt[i]);
    /* Accumulate to output buffer */
    ObitFArrayAdd (inImage->image, outImage->image, outImage->image);
 } /* end loop over input planes */

  /* Write output */
  plane[0] = outPlane;
  if ((ObitImagePutPlane (outImage, NULL, plane, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR writing image %s", 
		   routine, outImage->name);
    return;
  }
    
  /* Close - free buffers */
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
  if ((ObitImageClose (outImage, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, outImage->name);
    return;
  }
  /* Free image buffer if not memory resident */
  if (outImage->mySel->FileType!=OBIT_IO_MEM) 
    outImage->image = ObitFArrayUnref(outImage->image);
} /* end  ObitImageUtilInterp3 */

/**
 * Make antenna primary beam correction to an image based on the pointing
 * position in another image.
 * For frequencies < 1 GHz uses the VLA polynomial gain curves,
 * for higher frequencies, it uses a jinc function based on the antenna size.
 * \param inImage  Image to be corrected
 *     List items:
 * \li doTab OBIT_bool If True and a tabulated beam available, use it
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
  olong blc[IM_MAXDIM], trc[IM_MAXDIM];
  gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  odouble Freq;
  gboolean OK, isMeerKAT;
  ObitBeamShape *bs=NULL;
  ObitImageDesc *inDesc;
  olong nThreads, nTh, nrow, nrowPerThread, lorow, hirow;
  PBCorFuncArg **threadArgs;
  gchar *routine = "ObitImageUtilPBCorr";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(pntImage));
  g_assert (ObitImageIsA(outImage));

  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;

  if (antSize<0.01) antSize = 25.0; /* default antenna size */

  /* Beam shape */
  bs = ObitBeamShapeCreate ("BS", (ObitImage*)pntImage, 0.0, antSize, TRUE);
  /* Diagnostics */
  if (err->prtLv>=2) {
    if (bs->doTab) Obit_log_error(err, OBIT_InfoErr, "Using Tabulated Beam");
    if (bs->doVLITE) Obit_log_error(err, OBIT_InfoErr, "Using VLITE Beam");
    isMeerKAT = !strncmp(pntImage->myDesc->teles, "MeerKAT",7); /* MeerKAT? */
    if (isMeerKAT) Obit_log_error(err, OBIT_InfoErr, "Using MeerKAT Beam");
  }

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

  /* Set up - get frequency, default 1 GHz*/
  Freq = ObitImageMFGetPlaneFreq(inImage);
  if (Freq<=0.0) Freq = 1.0e9;
  ObitBeamShapeSetFreq(bs, Freq);  /* Set frequency */
  /* which beam model to use */
  inDesc = inImage->myDesc; /* Input descriptor */

  /* Initialize Threading */
  nThreads = MakePBCorFuncArgs (outImage->thread, bs, inDesc, 
				outImage->image, TRUE, 
				err, &threadArgs);

  /* Divide up work */
  nrow = outImage->myDesc->inaxes[1];
  nrowPerThread = nrow/nThreads;
  nTh = nThreads;
  if (nrow<64) {nrowPerThread = nrow; nTh = 1;}
  /* No fewer than 64 rows per thread */
  if ((nrowPerThread)<64) {
    nTh = (olong)(0.5 + nrow/64.0);
   if (nTh>0)  nrowPerThread = nrow/nTh;
  }
  if (nTh<=0) {nrowPerThread = nrow; nTh = 1;}
  lorow = 1;   hirow = nrowPerThread;  hirow = MIN (hirow, nrow);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hirow = nrow;  /* Make sure do all */
    threadArgs[i]->first   = lorow;
    threadArgs[i]->last    = hirow;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which row */
    lorow += nrowPerThread;
    hirow += nrowPerThread;
    hirow = MIN (hirow, nrow);
  }

  /* Do operation */
  OK = ObitThreadIterator (inImage->thread, nTh, 
			   (ObitThreadFunc)ThreadImagePBCor,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
  /* */

  /* Free local objects */
  KillPBCorFuncArgs(nThreads, threadArgs);

  /* Multiply by input - NB Beam inverted */
  ObitFArrayMul(outImage->image, inImage->image, outImage->image);

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
  /* Cleanup */
  bs = ObitBeamShapeUnref(bs);

} /* end  ObitImageUtilPBCorr */

/**
 * Multiply antenna primary beam pattern by an image based on the pointing
 * position in another image.
 * For frequencies < 1 GHz uses the VLA polynomial gain curves,
 * for higher frequencies, it uses a jinc function based on the antenna size.
 * \param inImage  Image to be corrected
 *     List items:
 * \li doTab OBIT_bool If True and a tabulated beam available, use it
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
  olong blc[IM_MAXDIM], trc[IM_MAXDIM];
  gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  odouble Freq;
  olong nThreads, nTh, nrow, nrowPerThread, lorow, hirow;
  gboolean  OK, isMeerKAT;
  ObitBeamShape *bs=NULL;
  ObitImageDesc *inDesc;
  PBCorFuncArg **threadArgs;
  gchar *routine = "ObitImageUtilPBApply";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(pntImage));
  g_assert (ObitImageIsA(outImage));

  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;

  if (antSize<0.01) antSize = 25.0; /* default antenna size */

  /* Beam shape */
  bs = ObitBeamShapeCreate ("BS", (ObitImage*)pntImage, 0.0, antSize, TRUE);
  /* Diagnostics */
  if (err->prtLv>=2) {
    if (bs->doTab) Obit_log_error(err, OBIT_InfoErr, "Using Tabulated Beam");
    if (bs->doVLITE) Obit_log_error(err, OBIT_InfoErr, "Using VLITE Beam");
    isMeerKAT = !strncmp(pntImage->myDesc->teles, "MeerKAT",7); /* MeerKAT? */
    if (isMeerKAT) Obit_log_error(err, OBIT_InfoErr, "Using MeerKAT Beam");
  }

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

  /* Set up - get frequency, default 1 GHz*/
  Freq = ObitImageMFGetPlaneFreq(inImage);
  if (Freq<=0.0) Freq = 1.0e9;
  ObitBeamShapeSetFreq(bs, Freq);  /* Set frequency */

  inDesc = inImage->myDesc; /* Input descriptor */

  /* Initialize Threading */
  nThreads = MakePBCorFuncArgs (outImage->thread, bs, inDesc, 
				outImage->image, FALSE, 
				err, &threadArgs);

  /* Divide up work */
  nrow = outImage->myDesc->inaxes[1];
  nrowPerThread = nrow/nThreads;
  nTh = nThreads;
  if (nrow<64) {nrowPerThread = nrow; nTh = 1;}
  /* No fewer than 64 rows per thread */
  if ((nrowPerThread)<64) {
    nTh = (olong)(0.5 + nrow/64.0);
   if (nTh>0)  nrowPerThread = nrow/nTh;
  }
  if (nTh<=0) {nrowPerThread = nrow; nTh = 1;}
  lorow = 1;   hirow = nrowPerThread;  hirow = MIN (hirow, nrow);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hirow = nrow;  /* Make sure do all */
    threadArgs[i]->first   = lorow;
    threadArgs[i]->last    = hirow;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which row */
    lorow += nrowPerThread;
    hirow += nrowPerThread;
    hirow = MIN (hirow, nrow);
  }

  /* Do operation */
  OK = ObitThreadIterator (inImage->thread, nTh, 
			   (ObitThreadFunc)ThreadImagePBCor,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
  /* */

  /* Free local objects */
  KillPBCorFuncArgs(nThreads, threadArgs);

  /* Multiply by input */
  ObitFArrayMul(outImage->image, inImage->image, outImage->image);

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

  /* Cleanup */
  bs = ObitBeamShapeUnref(bs);

} /* end  ObitImageUtilPBApply */

/**
 * Make an image of the antenna primary beam pattern based on the pointing
 * position in an image.  Honors existing blc,trc.
 * For frequencies < 1 GHz uses the VLA polynomial gain curves,
 * for higher frequencies, it uses a jinc function based on the antenna size.
 * \param pntImage Image with pointing position
 *     List items:
 * \li doTab    OBIT_bool If True and a tabulated beam available, use it
 * \li doInvert OBIT_bool If True invert beam gain (== correction)
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
  olong blc[IM_MAXDIM], trc[IM_MAXDIM];
  gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat ftemp, fblank = ObitMagicF();
  odouble Freq;
  olong nThreads, nTh, nrow, nrowPerThread, lorow, hirow;
  gboolean OK, isMeerKAT, doTab=FALSE, doInvert=FALSE;
  ObitInfoType type;
  ObitImageDesc *outDesc;
  ObitBeamShape *bs=NULL;
  PBCorFuncArg **threadArgs;
  gchar *routine = "ObitImageUtilPBImage";
 
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(pntImage));
  g_assert (ObitImageIsA(outImage));

  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
  dim[0] = 7;
  pntImage->extBuffer = FALSE;

  if (antSize<0.01) antSize = 25.0; /* default antenna size */
  ObitInfoListGetTest(pntImage->info, "doInvers", &type, dim, &doInvert); /* Invert? */
  ObitInfoListGetTest(pntImage->info, "doTab", &type, dim, &doTab); 

  /* Get pointing position */
  pntImage->extBuffer = TRUE;  /* Don't need buffer */
  if ((ObitImageOpen (pntImage, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, pntImage->name);
    return;
  }
 if ((ObitImageClose (pntImage, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, pntImage->name);
    return;
  }
  pntImage->extBuffer = FALSE;

  /* Beam shape */
  bs = ObitBeamShapeCreate ("BS", (ObitImage*)pntImage, minGain, antSize, TRUE);
  /* Diagnostics */
  if (err->prtLv>=2) {
    if (bs->doTab) Obit_log_error(err, OBIT_InfoErr, "Using Tabulated Beam");
    if (bs->doVLITE) Obit_log_error(err, OBIT_InfoErr, "Using VLITE Beam");
    isMeerKAT = !strncmp(pntImage->myDesc->teles, "MeerKAT",7); /* MeerKAT? */
    if (isMeerKAT) Obit_log_error(err, OBIT_InfoErr, "Using MeerKAT Beam");
    ObitErrLog(err); 
  }

  /* Do I/O by plane and all of plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListPut (outImage->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
  dim[0] = 7;
  for (i=0; i<5; i++) blc[i+2] = trc[i+2] = outPlane[i];
  ObitInfoListPut (outImage->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (outImage->info, "TRC", OBIT_long, dim, trc, err);
  outImage->myDesc->plane = outPlane[0];
  outImage->extBuffer = FALSE;


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

  /* Set up - get frequency, default 1 GHz*/
  Freq = ObitImageMFGetPlaneFreq(outImage);
  if (Freq<=0.0) Freq = 1.0e9;
  ObitBeamShapeSetFreq(bs, Freq);  /* Set frequency */
  /* Diagnostics */
  if (err->prtLv>=2) 
    Obit_log_error(err, OBIT_InfoErr, "Using Frequency %lg plane %d",
		   Freq,outPlane[0]);

  /* which beam model to use */
  outDesc = outImage->myDesc; /* output descriptor */

  /* Initialize Threading */
  nThreads = MakePBCorFuncArgs (outImage->thread, bs, outDesc, 
				outImage->image, doInvert, 
				err, &threadArgs);

  /* Divide up work */
  nrow = outImage->myDesc->inaxes[1];
  nrowPerThread = nrow/nThreads;
  nTh = nThreads;
  if (nrow<64) {nrowPerThread = nrow; nTh = 1;}
  /* No fewer than 64 rows per thread */
  if ((nrowPerThread)<64) {
    nTh = (olong)(0.5 + nrow/64.0);
   if (nTh>0)  nrowPerThread = nrow/nTh;
  }
  if (nTh<=0) {nrowPerThread = nrow; nTh = 1;}
  lorow = 1;   hirow = nrowPerThread;  hirow = MIN (hirow, nrow);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hirow = nrow;  /* Make sure do all */
    threadArgs[i]->first   = lorow;
    threadArgs[i]->last    = hirow;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which row */
    lorow += nrowPerThread;
    hirow += nrowPerThread;
    hirow = MIN (hirow, nrow);
  }

  /* Do operation */
  OK = ObitThreadIterator (outImage->thread, nTh, 
			   (ObitThreadFunc)ThreadImagePBCor,
			   (gpointer**)threadArgs);

  /* Check for problems */
  if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
  /* */

  /* Free local objects */
  KillPBCorFuncArgs(nThreads, threadArgs);

  /* Clip - check invert */
  if (doInvert && (minGain>1.0e-5)) {
    ftemp = 1.0 / minGain;
    ObitFArrayClip(outImage->image, -ftemp, ftemp, fblank);
  } else if (minGain>1.0e-5) {  /* No invert */
    ObitFArrayInClip(outImage->image, -minGain, minGain, fblank);
  }

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
  /* Cleanup */
  bs = ObitBeamShapeUnref(bs);

} /* end  ObitImageUtilPBImage */

/**
 * Make an image of the effective antenna primary beam pattern for On-The-Fly (OTF)
 * "Aussie mode: imaging based on the pointing position in an image.
 * For frequencies < 1 GHz uses the VLA polynomial gain curves,
 * for higher frequencies, it uses a jinc function based on the antenna size.
 * \param pntImage Image with pointing position
 *     List items:
 * \li doTab OBIT_bool If True and a tabulated beam available, use it
 * \param outImage Image to be written.  Must be previously instantiated.
 *     List items:
 * \li noff OBIT_long number of offsets
 * \li RAoff  OBIT_float (?,1,1) RA offsets (deg) from pntImage, not corrected by cos(dec)
 * \li Decoff OBIT_float (?,1,1) Dec offsets (deg) from pntImage
 * \param outPlane Desired plane in outImage on planes 3-5; ignored if memOnly
 * \param antSize  Antenna size
 * \param minGain  Min. allowed antenna gain, lower values are blanked
 * \param err      Error stack, returns if not empty.
 */
void 
ObitImageUtilOTFBeam (ObitImage *pntImage, ObitImage *outImage, 
		      olong *outPlane, ofloat antSize, ofloat minGain, ObitErr *err)
{
  ObitIOSize IOBy;
  olong blc[IM_MAXDIM], trc[IM_MAXDIM];
  gint32 i, dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong ix, iy, ip,  pos[2],noff;
  ollong indx;
  ofloat inPixel[2], *out, *norma, *RAoff, *Decoff, normc, fblank = ObitMagicF();
  odouble RAPnt, DecPnt, Freq, ra, dec, dist;
  odouble offRA, offDec;
  ofloat pbf, equinox;
  gboolean bad, isMeerKAT;
  ObitImageDesc *outDesc;
  ObitBeamShape *bs=NULL;
  ObitFArray *norm=NULL;
  gchar *routine = "ObitImageUtilOTFBeam";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(pntImage));
  g_assert (ObitImageIsA(outImage));
  /* Get offsets */
  noff = 1; 
  ObitInfoListGetTest(outImage->info, "noff", &type, dim, (gpointer)&noff);
  RAoff = NULL;
  ObitInfoListGetP(outImage->info, "RAoff", &type, dim, (gpointer*)&RAoff);
  Decoff = NULL;
  ObitInfoListGetP(outImage->info, "Decoff", &type, dim, (gpointer*)&Decoff);
  Obit_return_if_fail((RAoff && Decoff), err, "%s: RAoff or Decoff not given", routine);

  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;

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
  if ((ObitImageClose (pntImage, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s", 
		   routine, pntImage->name);
    return;
  }
  pntImage->extBuffer = FALSE;  /* May need buffer later */

   /* Beam shape */
  bs = ObitBeamShapeCreate ("BS", (ObitImage*)pntImage, minGain, antSize, TRUE);
  /* Diagnostics */
  if (err->prtLv>=2) {
    if (bs->doTab) Obit_log_error(err, OBIT_InfoErr, "Using Tabulated Beam");
    if (bs->doVLITE) Obit_log_error(err, OBIT_InfoErr, "Using VLITE Beam");
    isMeerKAT = !strncmp(pntImage->myDesc->teles, "MeerKAT",7); /* MeerKAT? */
    if (isMeerKAT) Obit_log_error(err, OBIT_InfoErr, "Using MeerKAT Beam");
  }

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

  /* Precess pointing position if needed */
  if ((abs(equinox-1950.0)<0.01) && 
      (abs(outImage->myDesc->equinox-2000.0)<0.01))
    ObitSkyGeomBtoJ (&RAPnt, &DecPnt);
  else if ((abs(equinox-2000.0)<0.01) && 
      (abs(outImage->myDesc->equinox-1950.0)<0.01))
    ObitSkyGeomJtoB (&RAPnt, &DecPnt);
  RAPnt  *= DG2RAD;
  DecPnt *= DG2RAD;


   /* Get output aray pointer */
  pos[0] = pos[1] = 0;
  out = ObitFArrayIndex (outImage->image, pos);
  /* Zero output image */
  ObitFArrayFill(outImage->image, 0.0);
  /* Normalization array */
  outDesc = outImage->myDesc; /* output descriptor */
  norm = ObitFArrayCreate("Norm", 2, outDesc->inaxes);
  ObitFArrayFill(norm, 0.0);
  norma = ObitFArrayIndex (norm, pos);

  /* Set up - get frequency, default 1 GHz*/
  outDesc->plane = outPlane[0];
  Freq = ObitImageMFGetPlaneFreq(outImage);
  if (Freq<=0.0) Freq = 1.0e9;
  ObitBeamShapeSetFreq(bs, Freq);  /* Set frequency */

 /* Loop over pointings */
  for (ip=0; ip<noff; ip++) {
    /* offset pointing */
    offRA  = RAPnt + RAoff[ip]*DG2RAD/cos(outDesc->crval[outDesc->jlocr]*DG2RAD);
    offDec = DecPnt + Decoff[ip]*DG2RAD;
    offRA *= RAD2DG; offDec *= RAD2DG;

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
			      offRA, offDec,
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
	dist = ObitBeamShapeAngle(bs, ra, dec, 0.0);
	
	/* primary beam correction */
	pbf = ObitBeamShapeGainSym (bs, dist);
	/* Clip by minGain */
	if (pbf<minGain) pbf = fblank;
      
	/* Save beam image */
	if (pbf!=fblank)  {out[indx] += pbf; norma[indx] += 1.0;}
	
      } /* end loop over x */
    } /* end loop over y */
  } /* end loop over offset */
  
    /* Normalize  */
  ObitFArrayDiv(outImage->image, norm, outImage->image);
  indx = outDesc->inaxes[0]*(1+outDesc->inaxes[1]/2) + (1+outDesc->inaxes[0]/2);
  normc = 1.0 / out[indx];
  ObitFArraySMul(outImage->image, normc);
  
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

  /* Cleanup */
  bs = ObitBeamShapeUnref(bs);

  } /* end  ObitImageUtilOTFBeam */

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
    maxfr = maxfr * sqrt(MaxBL / MaxW);
    /*  Assume circular clean beam - Undisturbed field size */
    *Radius =  maxfr*RAD2AS / (*Cells) + 0.5;
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
  olong pos[5]={0,0,0,0,0}, naxis[4] = {1,1,1,1};
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
  
  /* get axis dimensions */
  if (inArray->ndim>=1) naxis[0] = MAX (1, inArray->naxis[0]);
  if (inArray->ndim>=2) naxis[1] = MAX (1, inArray->naxis[1]);
  if (inArray->ndim>=3) naxis[2] = MAX (1, inArray->naxis[2]);
  if (inArray->ndim>=4) naxis[3] = MAX (1, inArray->naxis[3]);

  /* Loop over fourth dimension writing planes */
  for (j=0; j<naxis[3]; j++) {
    /* Loop over third */
    for (i=0; i<naxis[2]; i++) {
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

  /* Do I/O by plane (but keep any blc, trc) */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListPut (inImage->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
  inImage->extBuffer = FALSE;  /* Make sure it has buffer */

  /* Need statistics? */
  if (quant<=0.0) {

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
  /* Index as well */
  ObitImageDescIndex(out->myDesc);
  outDesc = out->myDesc;
  outDesc->bitpix = bitpix;
  outDesc->minval = minMin;
  outDesc->maxval = maxMax;
  g_snprintf (outDesc->origin, IMLEN_VALUE-1, "Generated by Obit");
  /* Save quantization value */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(outDesc->info, "Quant", OBIT_float, dim, &quant);

  /* Copy any descriptor info data */
  out->myDesc->info = ObitInfoListCopyData (inImage->myDesc->info, out->myDesc->info);

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
  outDesc->altCrpix = inDesc->altCrpix/incr;  /* Alt (vel) reference pixel */

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
  olong i, ip, np, pln[5];
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

  /* How many input planes */
  np = inDesc->inaxes[2];

  /* Read input plane */
  /* Set blc, trc */
  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<3; i++) trc[i] = inDesc->inaxes[i];
  for (i=3; i<IM_MAXDIM; i++) trc[i] = 1;
  
  dim[0] = 1;
  ObitInfoListAlwaysPut (in->info, "IOBy", OBIT_long, dim, &IOBy);
  dim[0] = 7;
  ObitInfoListAlwaysPut (in->info, "BLC", OBIT_long, dim, blc); 
  ObitInfoListAlwaysPut (in->info, "TRC", OBIT_long, dim, trc);


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
  pln[0] = pln[1] = pln[2] = pln[3] = pln[4] = 1;
  /* Loop over input planes */
  for (ip=0; ip<np; ip++) {
    pln[0] = ip+1;  pln[1] = pln[2] = pln[3] = pln[4] = 1;
    ObitImageGetPlane (in, NULL, pln, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    pln[0] = plane[0]+ip; 
    for (i=1; i<5; i++) pln[i] = plane[i];
    ObitImagePutPlane (out, in->image->array, pln, err);
    if (err->error) Obit_traceback_msg (err, routine, out->name);
  } /* end loop copying planes */
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
  olong inblc[IM_MAXDIM], blc[IM_MAXDIM], intrc[IM_MAXDIM], trc[IM_MAXDIM];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gchar *routine = "ObitImageUtilInsertPlane";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageIsA(in));
  g_assert (ObitImageIsA(out));
  g_assert (plane!=NULL);

  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
  for (i=0; i<IM_MAXDIM; i++) inblc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) intrc[i] = 0;

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
 * \param doGrid    If true force 2D images to common grid
 * \param nchavg    How many uv channels to average per image channel.
 *                  Ignored if uv data has multiple IFs.
 */
void 
ObitImageUtilUV2ImageDesc(ObitUVDesc *UVDesc, ObitImageDesc *imageDesc, 
			  gboolean doGrid, olong nchavg)
{
  olong i, iaxis, nif, nfreq, nch;
  odouble sum;
  gchar *st1, *st2;
  gchar *today=NULL;
  gchar *origin="Obit", *sinstr="-SIN";

  /* error checks */
  g_assert (ObitUVDescIsA(UVDesc));
  g_assert (ObitImageDescIsA(imageDesc));
  
  /* Be sure UV descriptor is indexed */
  ObitUVDescIndex(UVDesc);

  /* Creation date today */
  today = ObitToday();
  strncpy (imageDesc->date, today, IMLEN_VALUE-1);
  if (today) g_free(today);
 
  /* 2D/3D stuff */
  if (imageDesc->do3D) {
    /* Position from uv reference and shift */
    ObitSkyGeomXYShift (UVDesc->crval[UVDesc->jlocr], 
			UVDesc->crval[UVDesc->jlocd],
			UVDesc->xshift, UVDesc->yshift, 
			imageDesc->crota[1],
			&imageDesc->crval[0], &imageDesc->crval[1]);
    imageDesc->xshift = UVDesc->xshift;
    imageDesc->yshift = UVDesc->yshift;
    imageDesc->xPxOff = 0.0;
    imageDesc->yPxOff = 0.0;
  } else { /* 2D */
    ObitImageUtilTwoDShift (UVDesc, imageDesc, doGrid);
  } /* end 2D */

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
  for (i=0; i<4; i++)  {   /* Replace blank with -SIN */
    if (st2[i+4]==' ') st1[i+4] = sinstr[i];
    else               st1[i+4] = st2[i+4];
  }
  st1[9] = 0;

  /* Reference pixel for 3D */
  if (imageDesc->do3D) { /* 3D */
    imageDesc->crpix[iaxis] = 1.0 + imageDesc->inaxes[iaxis] / 2.0;
  } else { /* 2D 
    imageDesc->crpix[iaxis] += 1.0 + imageDesc->inaxes[iaxis] / 2.0;*/
  }

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
  for (i=0; i<4; i++)  {   /* Replace blank with -SIN */
    if (st2[i+4]==' ') st1[i+4] = sinstr[i];
    else               st1[i+4] = st2[i+4];
  }
  st1[9] = 0;

  /* Reference pixel for 3D */
  if (imageDesc->do3D) { /* 3D */
    imageDesc->crpix[iaxis] = 1.0 + imageDesc->inaxes[iaxis] / 2.0;
  } else { /* 2D 
    imageDesc->crpix[iaxis] += 1.0 + imageDesc->inaxes[iaxis] / 2.0;*/
  }

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
    nch = (olong)(0.999 + (ofloat)nfreq / nchavg); /* how many output channels? */
    imageDesc->inaxes[iaxis] = nch;
    imageDesc->cdelt[iaxis] = nchavg*UVDesc->cdelt[UVDesc->jlocf]; 
    imageDesc->crpix[iaxis] = 1.0 + (nchavg-1.0) / 2.0; /* reference pixel */
    /* average first nchavg uv channels out output reference frequency */
    sum = 0.0;
    for (i=0; i<nchavg; i++) sum += UVDesc->freqArr[i];
    imageDesc->crval[iaxis] = sum / (nchavg);
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
  strncpy (imageDesc->origin, origin,         IMLEN_VALUE-1);
  strncpy (imageDesc->bunit,  "JY/BEAM ",     IMLEN_VALUE-1);
  /* Set current date */
  ObitImageUtilCurDate (imageDesc->date, IMLEN_VALUE-1);

  /* Other */
  imageDesc->altCrpix     = UVDesc->altCrpix;
  imageDesc->altRef       = UVDesc->altRef;
  imageDesc->restFreq     = UVDesc->restFreq;
  imageDesc->VelReference = UVDesc->VelReference;
  imageDesc->VelDef       = UVDesc->VelDef;
  imageDesc->epoch        = UVDesc->epoch;
  imageDesc->equinox      = UVDesc->equinox;
  /* Pointing position */
  if ((fabs(UVDesc->obsra)>1.0e-5) || ((fabs(UVDesc->obsdec)>1.0e-5))) {
    imageDesc->obsra  = UVDesc->obsra;
    imageDesc->obsdec = UVDesc->obsdec;
  } else {
    imageDesc->obsra  = UVDesc->crval[UVDesc->jlocr];
    imageDesc->obsdec = UVDesc->crval[UVDesc->jlocd];
  }

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
 * \li "BLC"     OBIT_long (7) Lowest x,y,v pixel number selected [def 1,1,1]
 * \li "TRC"     OBIT_long (7) Highest x,y,v pixel number selected [def all]
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
  ObitIOCode   iretCode;
  ofloat       RMS, maxF, newVal, *Parms=NULL, fblank =  ObitMagicF();
  ofloat       vel, minAllow, FWHM;
  olong         iplane;
  ObitFArray   *accum=NULL, *accumWt=NULL;
  ObitFArray   *ConvFn=NULL;
  ObitImage    *scrImage=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM], blc0[IM_MAXDIM], trc[IM_MAXDIM], trc0[IM_MAXDIM];
  gchar        *today=NULL;
  olong        i, pos[IM_MAXDIM], Cen[2], temp, naxis[2];
  gboolean     odd;
  gchar        *routine = "ObitImageUtilVel";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));

  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
  for (i=0; i<IM_MAXDIM; i++) blc0[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc0[i] = 0;

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
  /* Not needed?ObitInfoListAlwaysPut (inImage->info, "TRC", OBIT_long, dim, trc);*/

  /* Blank or zero for out of range points? */
  if (Parms[2]<=0.5) newVal = fblank;
  else newVal = 0.0;

  /* Open input image */
  iretCode = ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Make scratch copy for convolution */
  scrImage = newObitImageScratch (inImage, err);
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
  ObitConvUtilConv (inImage, ConvFn, FALSE, FALSE, 1.0, scrImage, err);
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
  ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  /* Write plane */
  ObitImageWrite(outImage, accum->array, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  ObitImageClose (outImage, err);
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
 * \li "BLC"     OBIT_long (7) Lowest x,y,v pixel number selected [def 1,1,1]
 * \li "TRC"     OBIT_long (7) Highest x,y,v pixel number selected [def all]
 * \li "inc"     OBIT_long (7) Pixel increment on each axis [def all 1]
 * \param outImage Image to be written.  Must be previously instantiated.
 * \param err      Error stack, returns if not empty.
 */
void 
ObitImageUtilSelCopy (ObitImage *inImage, ObitImage *outImage, ObitErr *err)
{
  ObitIOCode   iretCode;
  ofloat       tmp;
  olong        itemp,  plane[IM_MAXDIM-2] = {0,1,1,1,1};
  olong        lblc[MAXFARRAYDIM], ltrc[MAXFARRAYDIM], linc[MAXFARRAYDIM];
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         inc[IM_MAXDIM];
  olong         blc[IM_MAXDIM], blc0[IM_MAXDIM];
  olong         trc[IM_MAXDIM], trc0[IM_MAXDIM];
  gchar        *today=NULL;
  olong        i;
  gboolean     want;
  gchar        *routine = "ObitImageUtilSelCopy";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));

  for (i=0; i<IM_MAXDIM; i++) blc[i]  = 1;
  for (i=0; i<IM_MAXDIM; i++) blc0[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i]  = 0;
  for (i=0; i<IM_MAXDIM; i++) trc0[i] = 0;
  for (i=0; i<IM_MAXDIM; i++) inc[i]  = 0;
  for (i=0; i<MAXFARRAYDIM; i++) linc[i] = lblc[i] = ltrc[i] = 1;

  /* Control parameters */
  ObitInfoListGetTest(inImage->info, "inc", &type, dim, inc); 
  for (i=0; i<IM_MAXDIM; i++) {
    if (inc[i]<=0) inc[i] = 1;
    linc[i] = MAX (1, inc[i]);
  }
  for (i=0; i<IM_MAXDIM; i++) blc[i] = MAX(1, blc0[i]);
  ObitInfoListGetTest(inImage->info, "BLC", &type, dim, blc);
  for (i=0; i<IM_MAXDIM; i++) trc[i] =  MAX(1, trc0[i]);
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
  ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Loop over planes until hitting EOF */
  while (iretCode!= OBIT_IO_EOF) {
    /* Which plane */
    plane[0]++;
    for (i=0; i<IM_MAXDIM-3; i++) {
      if (plane[i] > inImage->myDesc->inaxes[i+2]) {
	plane[i+1]++;
	plane[i] = 1;
	if (i>inImage->myDesc->naxis-3) break;
      }	else break;
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
    ObitImageWrite(outImage, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    
  } /* End loop over input image */

  /* Close files */
  iretCode = ObitImageClose (inImage, err);
  ObitImageClose (outImage, err);
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
  olong blc[IM_MAXDIM], trc[IM_MAXDIM];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOSize IOBy;
  gchar *today=NULL;
  ofloat Lambda, dx, dy, dist, xcenter, ycenter, val, radius2, arg;
  olong i, FFTdim[2], cen[2], ix, iy, even;
  ObitFArray *inFArray=NULL, *outFArray=NULL, *inFArrayCopy=NULL, *maskArray=NULL;
  ObitCArray *inCArray=NULL, *outCArray=NULL;
  ObitFFT *FFTfor=NULL, *FFTrev=NULL;
  gchar *routine = "ObitImageUtilUVFilter";
    
  /* error checks */
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));

  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;

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
  /*dx = fabs((FFTdim[0]*0.25*inImage->myDesc->cdelt[0]/57.296) / Lambda);
    dy = fabs((FFTdim[1]*0.25*inImage->myDesc->cdelt[1]/57.296) / Lambda);*/
  dx = 0.25*Lambda / fabs(FFTdim[0]*inImage->myDesc->cdelt[0]/57.296);
  dy = 0.25*Lambda / fabs(FFTdim[1]*inImage->myDesc->cdelt[1]/57.296);

  /* Form mask */
  xcenter = (ofloat)(FFTdim[0]/2); 
  ycenter = (ofloat)(FFTdim[1]/2);
  radius2 = radius*radius;
  for (iy=0; iy<FFTdim[1]; iy++) {
    for (ix=0; ix<FFTdim[0]; ix++) {
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
  /* Even or odd */
  if ((2*(inFArrayCopy->naxis[0]/2))==inFArrayCopy->naxis[0]) even = 1;
  else                                                        even = 0;
  trc[0] = cen[0] + even + inFArrayCopy->naxis[0] / 2;
  blc[1] = cen[1] - inFArrayCopy->naxis[1] / 2; 
  if ((2*(inFArrayCopy->naxis[1]/2))==inFArrayCopy->naxis[1]) even = 1;
  else                                                        even = 0;
  trc[1] = cen[1] + even + inFArrayCopy->naxis[1] / 2;
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

/**
 * Write the contents of an ObitFArray in a rudimentary FITS image
 * \param array    Array to write
 * \param FITSFile Name of FITS file
 * \param FITSdisk FITS disk number
 * \param desc     If nonNULL, use to derive the output header
 * \param err      Error stack, returns if not empty.
 * \return ObitImage* of output 
 */
ObitImage* ObitImageUtilFArray2FITS (ObitFArray *array, 
				     gchar *FITSFile, olong FITSdisk,
				     ObitImageDesc *desc, ObitErr *err)
{
  ObitImage  *outImage=NULL;
  olong      i, blc[IM_MAXDIM], trc[IM_MAXDIM];
  gchar      *routine = "ObitImageUtilFArray2FITS";

  /* error checks */
  if (err->error) return outImage;
  g_assert (ObitFArrayIsA(array));
  g_assert (FITSFile!=NULL);

  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;

  /* Create basic output Image Object */
  outImage = newObitImage("Output FITS Image");
  
  /* define object */
  ObitImageSetFITS (outImage, OBIT_IO_byPlane, FITSdisk, FITSFile, 
		    blc, trc, err);
  if (err->error) Obit_traceback_val (err, routine, FITSFile, outImage);

  /* Copy any descriptive material */
  if (desc) ObitImageDescCopyDesc (desc, outImage->myDesc, err);
  if (err->error) Obit_traceback_val (err, routine, FITSFile, outImage);

  /* Set Size */
  outImage->myDesc->bitpix    = -32;
  outImage->myDesc->naxis     = 2;
  outImage->myDesc->inaxes[0] = array->naxis[0];
  outImage->myDesc->inaxes[1] = array->naxis[1];

  /* Write */
  Obit_log_error(err, OBIT_InfoErr, "Write FITS image %d/ %s", FITSdisk,FITSFile);
  ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  ObitImageWrite(outImage, array->array, err);
  ObitImageClose (outImage, err);
  if (err->error) Obit_traceback_val (err, routine, FITSFile, outImage);
  
  return outImage;
} /* end ObitImageUtilFArray2FITS */


/**
 * Shift an image by a fixed amount.
 * Input image is shifted by an FFT/phase ramp/FFT method to another grid
 * shifted by shift pixels.
 * Processes a single plane.
 * Ouptut descriptor mostly unmodified - assumed coordinates are valid.
 * \param inImage  Image to be shifted.
 * \param outImage Image to be written.  Must be previously instantiated.
 * \param shift    Shift in pixels, in pixel space;
 *                 positive -> features will appear at higher pixel numbers.
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created ObitImage.
 */
void ObitImageUtilShift (ObitImage *inImage, ObitImage *outImage, ofloat *shift, 
			 ObitErr *err)
{
  ObitIOSize IOBy;
  olong blc[IM_MAXDIM], trc[IM_MAXDIM];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitFFT *FFTfor=NULL, *FFTrev=NULL;
  ObitFArray *inFArray=NULL, *outFArray=NULL;
  ObitCArray *inCArray=NULL, *outCArray=NULL;
  olong i, FFTdim[2], cen[2], ix, iy;
  ofloat xcenter, ycenter, xfact, yfact, phase, norm;
  gboolean doShift;
  gchar *today=NULL;
  gchar *routine = "ObitImageUtilShift";
 
  /* error checks */
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));

  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
  doShift = (shift[0]!=0.0) || (shift[1]!=0.0);  /* Need shift? */

  /* Do I/O by plane and all of plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListAlwaysPut (inImage->info, "IOBy", OBIT_long, dim, &IOBy);
  ObitInfoListAlwaysPut (outImage->info,"IOBy", OBIT_long, dim, &IOBy);
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (inImage->info, "BLC",  OBIT_long, dim, blc); 
  ObitInfoListAlwaysPut (inImage->info, "TRC",  OBIT_long, dim, trc);

  /* Open input  */
  if ((ObitImageOpen (inImage, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, inImage->name);
    return;
  }

  /* Read input plane */
  if ((ObitImageRead (inImage,NULL , err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR reading image %s", 
		   routine, inImage->name);
    return;
  }

  /* FFT size */
  FFTdim[0] = ObitFFTSuggestSize(inImage->myDesc->inaxes[0]);
  FFTdim[1] = ObitFFTSuggestSize(inImage->myDesc->inaxes[1]);

  /* Create float arrays for FFT size */
  inFArray  = ObitFArrayCreate("input",  2, FFTdim);
  outFArray = ObitFArrayCreate("output", 2, FFTdim);
    
  /* Pad input into work FArray */
  norm = 1.0 / (FFTdim[0]*FFTdim[1]); /* Normalization factor */
  ObitFArrayPad(inImage->image, inFArray, norm);
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
  /* Free image buffer if not memory resident */
  if (inImage->mySel->FileType!=OBIT_IO_MEM) 
    inImage->image = ObitFArrayUnref(inImage->image);

  /* Forward FFT */
  ObitFFTC2C (FFTfor, inCArray, outCArray);
    
  /* Calculate any phase ramp - use inCArray for work space */
  if (doShift) {
    xcenter = (ofloat)(FFTdim[0]/2); 
    ycenter = (ofloat)(FFTdim[1]/2);
    xfact = -2.0 * G_PI * shift[0] / FFTdim[0];
    yfact = -2.0 * G_PI * shift[1] / FFTdim[1];
    for (iy=0; iy<FFTdim[1]; iy++) {
      for (ix=0; ix<FFTdim[0]; ix++) {
	phase = (ix-xcenter)*xfact + (iy-ycenter)*yfact;
	inCArray->array[iy*FFTdim[0]*2+ix*2]   = cosf(phase);
	inCArray->array[iy*FFTdim[0]*2+ix*2+1] = sinf(phase);
      } /* end loop in x */
    } /* end loop in y */
  }

  /* Swaparoonie */
  ObitCArray2DCenterFull (inCArray);

  /* Apply phase ramp */
  if (doShift) ObitCArrayMul (outCArray, inCArray,outCArray);

  /* Create FFT for back FFT */
  FFTrev = newObitFFT("Forward FFT", OBIT_FFT_Reverse, OBIT_FFT_FullComplex, 2, FFTdim);

  /* Back FFT */
  ObitFFTC2C(FFTrev, outCArray, inCArray);
    
  /* Extract Real */
  ObitCArrayReal (inCArray, outFArray);
  /* and God said "The center of an FFT will be at the corners" */
  ObitFArray2DCenter(outFArray);
    
  /* Open output */
  if ((ObitImageOpen (outImage, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, outImage->name);
    return;
  }

  /* Extract output portion */
  cen[0] = FFTdim[0]/2; cen[1] = FFTdim[1]/2;
  blc[0] = cen[0] - inFArray->naxis[0] / 2; 
  trc[0] = cen[0] - 1 + inFArray->naxis[0] / 2;
  blc[1] = cen[1] - inFArray->naxis[1] / 2; 
  trc[1] = cen[1] - 1 + inFArray->naxis[1] / 2;
  ObitImageUnref(outImage->image);
  outImage->image = ObitFArraySubArr(outFArray, blc, trc, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Blank output where input blanked */
  ObitFArrayBlank (outImage->image, inFArray, outImage->image);

  /* Copy descriptor */
  /*  ObitImageDescCopyDesc (inImage->myDesc, outImage->myDesc, err);*/
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  
  /* Creation date today */
  today = ObitToday();
  strncpy (outImage->myDesc->date, today, IMLEN_VALUE-1);
  if (today) g_free(today);

  /* Update descriptor */
  outImage->myDesc->maxval = -1.0e20;
  outImage->myDesc->minval =  1.0e20;
  outImage->myDesc->areBlanks = FALSE;
  /*  outImage->myDesc->crpix[0] += shift[0];*/
  /*  outImage->myDesc->crpix[1] += shift[1];*/
  if (!outImage->myDesc->do3D) {
    /*   outImage->myDesc->xPxOff -= shift[0];*/
    /*   outImage->myDesc->yPxOff -= shift[1];*/
  }

  /* Write */
  ObitImageWrite (outImage, outImage->image->array, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  ObitImageClose (outImage,err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  /* Free image buffer if not memory resident */
  if (outImage->mySel->FileType!=OBIT_IO_MEM) 
    outImage->image = ObitFArrayUnref(outImage->image);

  /* Cleanup */
  if (inImage->image)  inImage->image  = ObitFArrayUnref(inImage->image);
  if (outImage->image) outImage->image = ObitFArrayUnref(outImage->image);
  if (inFArray)     ObitFArrayUnref(inFArray);
  if (outFArray)    ObitFArrayUnref(outFArray);
  if (inCArray)     ObitCArrayUnref(inCArray);
  if (outCArray)    ObitCArrayUnref(outCArray);
  if (FFTfor)       ObitFFTUnref(FFTfor);
  if (FFTrev)       ObitFFTUnref(FFTrev);
}  /* end ObitImageUtilShift */

/**
 * Shift an MF image by a fixed amount.
 * Input image is shifted by an FFT/phase ramp/FFT method to another grid
 * shifted by shift pixels.
 * Loops over coarse planes of the MF Image.
 * \param inImage  Image to be shifted.
 * \param outImage Image to be written.  Must be previously instantiated.
 * \param shift    Shift in pixels, in pixel space;
 *                 positive -> features will appear at higher pixel numbers.
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created ObitImage.
 */
void ObitImageUtilMFShift (ObitImage *inImage, ObitImage *outImage, ofloat *shift, 
			   ObitErr *err)
{
  ObitIOSize IOBy;
  ObitIOCode retCode;
  olong blc[IM_MAXDIM], trc[IM_MAXDIM], plane[5] = {1,1,1,1,1};
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitFFT *FFTfor=NULL, *FFTrev=NULL;
  ObitFArray *inFArray=NULL, *outFArray=NULL;
  ObitCArray *inCArray=NULL, *outCArray=NULL;
  olong i, FFTdim[2], cen[2], ix, iy, nSpec, nTerm, iplane;
  ofloat xcenter, ycenter, xfact, yfact, phase, norm;
  gchar *today=NULL;
  gboolean doShift;
  gchar *routine = "ObitImageUtilMFShift";
 
  /* error checks */
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));
  Obit_return_if_fail(ObitImageMFIsA((ObitImageMF*)inImage), err,
		      "%s: input image %s not MF", 
		      routine, inImage->name);
  Obit_return_if_fail(ObitImageMFIsA((ObitImageMF*)outImage), err,
		      "%s: output image %s not MF", 
		      routine, outImage->name);
  
  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;

  /* MF Info */
  nSpec = ((ObitImageMF*)inImage)->nSpec;
  nTerm = ((ObitImageMF*)inImage)->maxOrder + 1;
  doShift = (shift[0]!=0.0) || (shift[1]!=0.0);  /* Need shift? */

  /* Do I/O by plane and all of plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListAlwaysPut (inImage->info, "IOBy", OBIT_long, dim, &IOBy);
  ObitInfoListAlwaysPut (outImage->info,"IOBy", OBIT_long, dim, &IOBy);
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (inImage->info,  "BLC",  OBIT_long, dim, blc); 
  ObitInfoListAlwaysPut (inImage->info,  "TRC",  OBIT_long, dim, trc);
  ObitInfoListAlwaysPut (outImage->info, "BLC",  OBIT_long, dim, blc); 
  ObitInfoListAlwaysPut (outImage->info, "TRC",  OBIT_long, dim, trc);
  
  /* Open input  */
  retCode = ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, inImage->name);
    return;
  }
  
  /* Open output */
  retCode = ObitImageOpen (outImage, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, outImage->name);
    return;
  }
  
  /* FFT size */
  FFTdim[0] = ObitFFTSuggestSize(inImage->myDesc->inaxes[0]);
  FFTdim[1] = ObitFFTSuggestSize(inImage->myDesc->inaxes[1]);
  
  /* Create FFT for full complex FFT */
  FFTfor = newObitFFT("Forward FFT", OBIT_FFT_Forward, OBIT_FFT_FullComplex, 2, FFTdim);
  
  /* Create FFT for back FFT */
  FFTrev = newObitFFT("Forward FFT", OBIT_FFT_Reverse, OBIT_FFT_FullComplex, 2, FFTdim);
  
  /* Create float arrays for FFT size */
  inFArray  = ObitFArrayCreate("input",  2, FFTdim);
  outFArray = ObitFArrayCreate("output", 2, FFTdim);
  
  /* Create complex arrays for FFT size */
  inCArray  = ObitCArrayCreate("input", 2, FFTdim);
  outCArray = ObitCArrayCreate("output", 2, FFTdim);
  
  /* Copy descriptor */
  ObitImageDescCopyDesc (inImage->myDesc, outImage->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  
  /* Creation date today */
  today = ObitToday();
  strncpy (outImage->myDesc->date, today, IMLEN_VALUE-1);
  if (today) g_free(today);
  
  /* Update descriptor */
  outImage->myDesc->maxval = -1.0e20;
  outImage->myDesc->minval =  1.0e20;
  outImage->myDesc->areBlanks = FALSE;
  outImage->myDesc->crpix[0] += shift[0];
  outImage->myDesc->crpix[1] += shift[1];
  if (!outImage->myDesc->do3D) {
    outImage->myDesc->xPxOff -= shift[0];
    outImage->myDesc->yPxOff -= shift[1];
  }

  /* Loop over coarse planes */
  for (iplane = 0; iplane<nSpec; iplane++) {
    plane[0] = 1 + nTerm + iplane;
    
    /* Read input plane */
    retCode = ObitImageGetPlane (inImage, NULL,  plane,err);
    if ((retCode != OBIT_IO_OK) || (err->error>0)) { /* error test */
      Obit_log_error(err, OBIT_Error, "%s: ERROR reading image %s", 
		     routine, inImage->name);
      return;
    }
    
    /* Pad input into work FArray */
    norm = 1.0 / (FFTdim[0]*FFTdim[1]); /* Normalization factor */
    ObitFArrayPad(inImage->image, inFArray, norm);
    /* and God said "The center of an FFT will be at the corners" */
    ObitFArray2DCenter (inFArray);
    /* Zero output FArray and use as imaginary part */
    ObitFArrayFill(outFArray, 0.0);
    /* Copy input to scratch CArray */
    ObitCArrayComplex(inFArray, outFArray, inCArray);
    
    /* Forward FFT */
    ObitFFTC2C (FFTfor, inCArray, outCArray);
    
    /* Calculate any phase ramp - use inCArray for work space */
    if (doShift) {
      xcenter = (ofloat)(FFTdim[0]/2); 
      ycenter = (ofloat)(FFTdim[1]/2);
      xfact = -2.0 * G_PI * shift[0] / FFTdim[0];
      yfact = -2.0 * G_PI * shift[1] / FFTdim[1];
      for (iy=0; iy<FFTdim[1]; iy++) {
	for (ix=0; ix<FFTdim[0]; ix++) {
	  phase = (ix-xcenter)*xfact + (iy-ycenter)*yfact;
	  inCArray->array[iy*FFTdim[0]*2+ix*2]   = cosf(phase);
	  inCArray->array[iy*FFTdim[0]*2+ix*2+1] = sinf(phase);
	} /* end loop in x */
      } /* end loop in y */
    }
    /* Swaparoonie */
    ObitCArray2DCenterFull (inCArray);
    
    /* Apply any phase ramp */
    if (doShift) ObitCArrayMul (outCArray, inCArray, outCArray);
    
    /* Back FFT */
    ObitFFTC2C(FFTrev, outCArray, inCArray);
    
    /* Extract Real */
    ObitCArrayReal (inCArray, outFArray);
    /* and God said "The center of an FFT will be at the corners" */
    ObitFArray2DCenter(outFArray);
    
    /* Extract output portion */
    cen[0] = FFTdim[0]/2; cen[1] = FFTdim[1]/2;
    blc[0] = cen[0] - inFArray->naxis[0] / 2; 
    trc[0] = cen[0] - 1 + inFArray->naxis[0] / 2;
    blc[1] = cen[1] - inFArray->naxis[1] / 2; 
    trc[1] = cen[1] - 1 + inFArray->naxis[1] / 2;
    ObitImageUnref(outImage->image);
    outImage->image = ObitFArraySubArr(outFArray, blc, trc, err);
    if (err->error) Obit_traceback_msg (err, routine, inImage->name);
    
    /* Blank output where input blanked */
    ObitFArrayBlank (outImage->image, inFArray, outImage->image);
    
    /* Write */
    retCode = ObitImagePutPlane (outImage, outImage->image->array, plane, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  } /* end loop over planes */

  /* Close files */
  ObitImageClose (inImage, err) ;
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  /* Free image buffer if not memory resident */
  if (inImage->mySel->FileType!=OBIT_IO_MEM) 
    inImage->image = ObitFArrayUnref(inImage->image);
  
  ObitImageClose (outImage,err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  /* Free image buffer if not memory resident */
  if (outImage->mySel->FileType!=OBIT_IO_MEM) 
    outImage->image = ObitFArrayUnref(outImage->image);

  /* Form combined image */
  ObitImageMFCombine ((ObitImageMF*)outImage, doShift, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Cleanup */
  if (inFArray)     ObitFArrayUnref(inFArray);
  if (outFArray)    ObitFArrayUnref(outFArray);
  if (inCArray)     ObitCArrayUnref(inCArray);
  if (outCArray)    ObitCArrayUnref(outCArray);
  if (FFTfor)       ObitFFTUnref(FFTfor);
  if (FFTrev)       ObitFFTUnref(FFTrev);
}  /* end ObitImageUtilMFShift */

olong ObitImageUtilBufSize (ObitUV *inUV)
{
  olong NPIO, nThreads;
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(inUV->thread));
  NPIO = nThreads * 
    (olong) (0.5 + MAX (2.0, 2.0/(inUV->myDesc->lrec*4.0/(1024.0*1024.0))));
  return NPIO;
} /* end ObitImageUtilBufSize */


/**
 * Fill in an image Descriptor shift information for 2D imaging
 * Shift is the closest to the shift in UVDesc such that pixels are 
 * on the same grid as the main tangent plane image.
 * \param UVDesc    Input UV Descriptor.
 * \param imageDesc Output image Descriptor
 * \param onGrid    If TRUE force alignment with a regular grid.
 */
void ObitImageUtilTwoDShift (ObitUVDesc *UVDesc, ObitImageDesc *imageDesc,
			     gboolean onGrid)
{
  ofloat  xpix, ypix, xcrpix, ycrpix, xshift, yshift;
  odouble ra, dec, rac, decc;

  /* First calculate the reference pixel using the shift in UVDesc */
  ObitSkyGeomXYShift (UVDesc->crval[UVDesc->jlocr], 
		      UVDesc->crval[UVDesc->jlocd],
		      UVDesc->xshift, UVDesc->yshift, imageDesc->crota[1],
		      &ra, &dec);
  ObitSkyGeomShiftCRP (&UVDesc->ptype[UVDesc->ilocu][4], 
		       UVDesc->crval[UVDesc->jlocr], UVDesc->crval[UVDesc->jlocd],
		       imageDesc->crota[1], ra, dec,
		       &xshift, &yshift);
  xcrpix = -xshift / imageDesc->cdelt[0];
  ycrpix = -yshift / imageDesc->cdelt[1];

  xcrpix += 1.0 + imageDesc->inaxes[0] / 2.0; /* offset to ref pix.*/
  ycrpix += 1.0 + imageDesc->inaxes[1] / 2.0;
  if (onGrid) {
    /* Get closest integral pixel */
    if (xcrpix>0.0)  
      xpix = (ofloat)((olong)(xcrpix+0.5));
    else
      xpix = (ofloat)((olong)(xcrpix-0.5));
    if (ycrpix>0.0)  
      ypix = (ofloat)((olong)(ycrpix+0.5));
    else
      ypix = (ofloat)((olong)(ycrpix-0.5));
  } else {  /* Not aligned with grid */
     xpix = xcrpix;
     ypix = ycrpix;
  }

  /* Get coordinate (rac, decc) of image at (xpix,ypix) */
  ObitSkyGeomWorldPos (xpix, ypix, ra, dec, xcrpix, ycrpix,
		       imageDesc->cdelt[0], imageDesc->cdelt[1], 
		       imageDesc->crota[1], &imageDesc->ctype[0][4],
		       &rac, &decc);

  /* Get shift to position (rac,decc) */
  ObitSkyGeomShiftXY (UVDesc->crval[UVDesc->jlocr], UVDesc->crval[UVDesc->jlocd],
		      imageDesc->crota[1], rac, decc,
		      &xshift, &yshift);
  imageDesc->xshift = xshift;
  imageDesc->yshift = yshift;

  /* Position  - as offset from UV (tangent) position */
  imageDesc->crval[0] = UVDesc->crval[UVDesc->jlocr];
  imageDesc->crval[1] = UVDesc->crval[UVDesc->jlocd];
  imageDesc->crpix[0] = xpix;
  imageDesc->crpix[1] = ypix;

  /* Correct version */
  imageDesc->xPxOff = xpix - (1.0 + imageDesc->inaxes[0] / 2.0);
  imageDesc->yPxOff = ypix - (1.0 + imageDesc->inaxes[1] / 2.0);

} /* end ObitImageUtilTwoDShift */

/**
 * Convert an ObitImage(MF) (TSpec CCs) to an ObitImageWB (Spec CCs)
 * Copies spectral planes and converts specified CC table
 * Tabulated spectrum fitted with spectrum weighting by primary beam
 * \param inImage  Input Image.
 *                 possible control parameter:
 * \li "dropNeg" OBIT_boolean if True, drop negative components [True] 
 * \param outImage Image to be written.  Must be previously defined.
 *                 Returned as ObitImageWB.
 * \param nTerm    Number of output Spectral terms, 2=SI, 3=also curve.
 * \param inCCVer  Input CCTable to convert, 0=> highest
 * \param outCCVer Output CCTable, 0=>1
 * \param startCC  First 1-rel component to convert
 * \param endCC    Last 1-rel component to convert, 0=> all
 * \param err      Error stack, returns if not empty.
 */
void ObitImageUtilT2Spec  (ObitImage *inImage, ObitImage **outImage, 
			   olong nTerm, olong *inCCVer, olong *outCCVer,
			   olong startCC, olong endCC, ObitErr *err)
{
  ObitImageWB *outImageWB=NULL;
  ObitImageDesc *imDesc=NULL;
  ObitIOSize IOBy;
  ObitIOCode retCode;
  ObitHistory *inHist=NULL, *outHist=NULL;
  olong blc[IM_MAXDIM], trc[IM_MAXDIM];
  olong inTerm, nSpec, norder;
  olong i, iplane, planeNo[5] = {1,1,1,1,1};
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *today=NULL;
  gchar *routine = "ObitImageUtilT2Spec";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageIsA(inImage));

  /* Do I/O by plane and all of plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListAlwaysPut (inImage->info,  "IOBy", OBIT_long, dim, &IOBy);

  /* Select portion of image */
  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;

  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (inImage->info,  "BLC",  OBIT_long, dim, blc); 
  ObitInfoListAlwaysPut (inImage->info,  "TRC",  OBIT_long, dim, trc);

  /* Size of spectra, number spectral terms */
  nSpec = 1;
  ObitInfoListGetTest(inImage->myDesc->info, "NSPEC", &type, dim, &nSpec);
  inTerm = 1;
  ObitInfoListGetTest(inImage->myDesc->info, "NTERM", &type, dim, &inTerm);
  
  /* Open input  */
  retCode = ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, inImage->name);
    return;
  }

  /* Make sure this is an ObitImageMF - freq axis = "SPECLNMF" */
  imDesc = inImage->myDesc;
  Obit_return_if_fail((!strncmp (imDesc->ctype[imDesc->jlocf],"SPECLNMF", 8)), err, 
		      "%s: Image %s NOT an ObitImageMF - no SPECLNMF axis", 
		      routine, inImage->name);

  /* Generate output descriptor */
  ObitImageDescCopyDesc (imDesc, (*outImage)->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  
  /* Creation date today */
  today = ObitToday();
  strncpy ((*outImage)->myDesc->date, today, IMLEN_VALUE-1);
  if (today) g_free(today);
  
  /* Update descriptor */
  (*outImage)->myDesc->bitpix = -32;
  (*outImage)->myDesc->maxval = -1.0e20;
  (*outImage)->myDesc->minval =  1.0e20;
  (*outImage)->myDesc->areBlanks = FALSE;
  (*outImage)->myDesc->naxis  = imDesc->naxis;
  for (i=0; i<IM_MAXDIM; i++) (*outImage)->myDesc->inaxes[i] = 1;
  (*outImage)->myDesc->inaxes[0] = imDesc->inaxes[0];
  (*outImage)->myDesc->inaxes[1] = imDesc->inaxes[1];
  (*outImage)->myDesc->inaxes[2] = 1;
  /* Convert to ObitImageWB */
  outImageWB = ObitImageWBFromImage (*outImage, nTerm-1, err);
  if (err->error) Obit_traceback_msg (err, routine, (*outImage)->name);
  /* Do I/O by plane and all of plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListAlwaysPut (outImageWB->info, "IOBy", OBIT_long, dim, &IOBy);
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (outImageWB->info, "BLC",  OBIT_long, dim, blc); 
  ObitInfoListAlwaysPut (outImageWB->info, "TRC",  OBIT_long, dim, trc);
  (*outImage) = ObitImageUnref(*outImage);
  (*outImage) = (ObitImage*)ObitImageWBRef(outImageWB);
  
  /* use same data buffer on input and output 
     so don't assign buffer for output */
  (*outImage)->extBuffer = TRUE;

  /* Open output */
  retCode = ObitImageOpen ((*outImage), OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, (*outImage)->name);
    return;
  }

  /* Copy first two planes */
  for (iplane=0; iplane<2; iplane++) {
    /* Write it */
    planeNo[0] = iplane+1; 
    ObitImageGetPlane (inImage, NULL, planeNo, err);
    ObitImagePutPlane ((*outImage), inImage->image->array, 
		       planeNo, err);
    if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  } /* end loop writing planes */

  /* Copy/convert table */
  ObitTableCCUtilT2Spec (inImage, (ObitImageWB*)(*outImage), nTerm, inCCVer, outCCVer, 
			 startCC, endCC, err);
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);

  /* Copy any history  unless Scratch */
  if (!inImage->isScratch && !(*outImage)->isScratch) {
    inHist  = newObitDataHistory((ObitData*)inImage, OBIT_IO_ReadOnly, err);
    outHist = newObitDataHistory((ObitData*)(*outImage), OBIT_IO_WriteOnly, err);
    outHist = ObitHistoryCopy (inHist, outHist, err);
    if (err->error) Obit_traceback_msg (err, routine, inImage->name);
    inHist  = ObitHistoryUnref(inHist);
    outHist = ObitHistoryUnref(outHist);
  }
  /* Make sure output closed */
  ObitImageClose ((*outImage),err);
  if (err->error) Obit_traceback_msg (err, routine, (*outImage)->name);
  (*outImage)->extBuffer = FALSE;  /* May need buffer later */

  /* Open output */
  retCode = ObitImageOpen ((*outImage), OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, (*outImage)->name);
    return;
  }
  /* Relabel output frequency axis NO!
     strncpy((*outImage)->myDesc->ctype[(*outImage)->myDesc->jlocf], "FREQ    ", 8); */

  /* Cleanup output Descriptor list */
  (*outImage)->myDesc->info = ObitInfoListUnref((*outImage)->myDesc->info);
  (*outImage)->myDesc->info = newObitInfoList();

  /* Add order to descriptor infolist */
  dim[0] = dim[1] = dim[2] = 1;
  norder = nTerm-1;
  ObitInfoListAlwaysPut((*outImage)->myDesc->info, "NTERM", OBIT_long, dim, &norder);
  /* Force update */
  (*outImage)->myStatus = OBIT_Modified;

  /* Close files */
  ObitImageClose (inImage, err) ;
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  /* Free image buffer if not memory resident */
  if (inImage->mySel->FileType!=OBIT_IO_MEM) 
    inImage->image = ObitFArrayUnref(inImage->image);

  ObitImageClose ((*outImage),err);
  if (err->error) Obit_traceback_msg (err, routine, (*outImage)->name);
  (*outImage)->extBuffer = FALSE;  /* May need buffer later */

  } /* end ObitImageUtilT2Spec */

/**
 * Fit Gaussian to a Beam.
 * Adapted from AIPS BMSHP.FOR, FITBM.FOR
 * \param beam  Beam image to fit.  
 *        Plane specified in info list member "PLANE" [def 1]
 * \param err   Error stack, returns if not empty.
 */
void ObitImageUtilFitBeam (ObitImage *beam, ObitErr *err)
{
  olong blc[2], trc[2], center[2], prtLv, iplane, plane[] = {1,1,1,1,1};
  ofloat peak, cellx, celly, fblank =  ObitMagicF();
  ofloat bmaj, bmin, bpa, dx, dy;
  olong itemp, nmodel, nparm, corner[2], fdim[2], hwid;
  ofloat Peak, RMSResid, peakResid, fluxResid, DeltaX,  DeltaY, parms[3], BeamTapr;
  ObitFArray *beamData = NULL, *beamCenter = NULL, *matx = NULL;
  ObitImageDesc *desc=(ObitImageDesc*)beam->myIO->myDesc;
  ObitImageFit *imFit=NULL;
  ObitFitRegion *reg=NULL;
  ObitFitModel **models = {NULL};
  ObitInfoType type;
  gboolean True=TRUE;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "ObitImageUtilFitBeam";

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageIsA(beam));

  /* plane number - first one with data */
  iplane = 1;
  if (desc->inaxes[2]>1)
    ObitInfoListGetTest(beam->info, "PLANE", &type, dim, &iplane);
  while (iplane<=desc->inaxes[2]) {
    plane[0] = iplane++;
    
    /* Read beam until non blanked plane */
    ObitImageGetPlane (beam, NULL, plane, err);
    if (err->error) goto cleanup;
    if (beam->image->array[10]!=fblank) break;  /* This one OK? */
  }
  /* test if have beam */
  Obit_return_if_fail((beam->image!=NULL), err, 
		      "%s: Problem with beam",  routine);
  /* Trim edges of beam */
  blc[0] = beam->image->naxis[0]/4;
  blc[1] = beam->image->naxis[1]/4;
  trc[0] = 3*beam->image->naxis[0]/4;
  trc[1] = 3*beam->image->naxis[1]/4;
  beamCenter = ObitFArraySubArr (beam->image, blc, trc, err);
  if (err->error) goto cleanup;

  /* Find center */
  peak = ObitFArrayMax (beamCenter, center);
  if (peak<0.1) goto cleanup;

  /* Close Beam */
  if ((ObitImageClose (beam, err)  
       != OBIT_IO_OK) || (err->error>0)) {  /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR closing image %s",  
		   routine, beam->name); 
    goto cleanup; 
  }

  /* Center better be 1.0 */
  if (fabs (peak-1.0) > 0.001) { 
    Obit_log_error(err, OBIT_Error,  
      		   "Beam peak (%f) not 1.0 for %s",  
      		   peak, beam->name); 
    goto cleanup; 
  } 

  /* Image info */
  cellx = beam->myDesc->cdelt[0];
  celly = beam->myDesc->cdelt[1];

  /* Fit Gaussian */
  /* Initial model */
  Peak   = 1.0;
  nparm  = 3;
  parms[0]  = parms[1] = 6.0; parms[2] = 0.0;
  nmodel = 1;
  /* Size of fitting region depends on BeamTapr, if any */
  hwid = 25;
  if (ObitInfoListGetTest(beam->myDesc->info, "BeamTapr", &type, dim, &BeamTapr)) {
    itemp = (olong)(0.5+BeamTapr/celly);
    hwid = MAX (hwid, itemp);
  }
  DeltaX = hwid+1; DeltaY = hwid+1;

  models = g_malloc0(sizeof(ObitFitModel*));
  models[0] = ObitFitModelCreate ("model", OBIT_FitModel_GaussMod, 
				  Peak, DeltaX, DeltaY, nparm, parms);
  /* Fitting region */
  corner[0] = center[0]+blc[0]-hwid; corner[1] = center[1]+blc[1]-hwid;
  fdim[0]   = fdim[1] = 2*hwid+1;
  peakResid = fluxResid = RMSResid = 0.0;
  reg    = ObitFitRegionCreate ("reg", corner, fdim, Peak, 
				RMSResid, peakResid, fluxResid, 
				nmodel, models);

  /* Fit in pixels/deg */
  imFit  = ObitImageFitCreate ("Fitter");
  /* Fix Peak, pos */
  dim[0] = dim[1] = dim[2] = dim[3] = 1;
  ObitInfoListAlwaysPut(imFit->info, "FixFlux",  OBIT_bool, dim, &True);
  ObitInfoListAlwaysPut(imFit->info, "FixPos",  OBIT_bool, dim, &True);

  /* Turn off fitting messages */
  prtLv = err->prtLv;
  err->prtLv = 0;
  ObitImageFitFit (imFit, beam, reg, err);
  err->prtLv = prtLv;
  if (err->error) goto cleanup;

  /* If different cells spacings */
  if (fabs(fabs(cellx)-fabs(celly)) > (0.01*fabs(celly))) {
    dx = reg->models[0]->parms[0]*fabs(cellx)*sin(reg->models[0]->parms[2]);
    dy =reg-> models[0]->parms[0]*fabs(celly)*cos(reg->models[0]->parms[2]);
    reg->models[0]->parms[0] = sqrt (dx*dx + dy*dy);
    reg->models[0]->parms[2] = atan2 (dx, dy);
    dx = reg->models[0]->parms[1]*fabs(cellx)*sin(reg->models[0]->parms[2]+G_PI/2);
    dy = reg->models[0]->parms[1]*fabs(celly)*cos(reg->models[0]->parms[2]+G_PI/2);
    reg->models[0]->parms[1] = sqrt (dx*dx + dy*dy);
  }

  bmaj = reg->models[0]->parms[0]*fabs(cellx); /* from pixels to deg */
  bmin = reg->models[0]->parms[1]*fabs(celly);
  bpa  = reg->models[0]->parms[2]*RAD2DG;

  /* add map rotation, force to +/- 90 deg. */
  bpa = bpa - beam->myDesc->crota[1];
  if (bpa > 90.0) bpa = bpa - 180.0;
  if (bpa < -90.0) bpa = bpa + 180.0;
  
  /* Give informative message about fitted size */
  Obit_log_error(err, OBIT_InfoErr, 
		 "Fitted beam %f x %f asec, PA = %f for %s plane %d", 
		 bmaj*3600.0, bmin*3600.0, bpa, beam->name, iplane);

  /* Save fitted values */
  beam->myDesc->beamMaj = bmaj;
  beam->myDesc->beamMin = bmin;
  beam->myDesc->beamPA  = bpa;

  /* Cleanup */
 cleanup: matx = ObitFArrayUnref (matx);
  beamData   = ObitFArrayUnref(beamData);
  beamCenter = ObitFArrayUnref(beamCenter);
  reg        = ObitFitRegionUnref(reg);
  imFit      = ObitImageFitUnref(imFit);
  if (models) {
    if (models[0]) models[0]  = ObitFitModelUnref(models[0]);
    g_free(models);
  }

  /* Free image buffer */
  beam->image = ObitFArrayUnref(beam->image);
  if (err->error) Obit_traceback_msg (err, routine, beam->name);

} /* end ObitImageUtilFitBeam */

/**
 * Fill an image with blank pixels (3rd dim only)
 * \param in        Input image to blank
 * \param err       Error stack, returns if not empty.
 */
void ObitImageUtilBlankFill (ObitImage* in, ObitErr* err)
{
  olong i, iplane, nplane;
  olong blc[IM_MAXDIM], trc[IM_MAXDIM];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat fblank = ObitMagicF();
  ObitIOSize IOBy = OBIT_IO_byPlane;
  gchar *routine="ObitImageUtilBlankFill";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageIsA(in));

  /* Access images a plane at a time */
  dim[0] = 1;
  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
  ObitInfoListPut (in->info, "IOBy", OBIT_long, dim, &IOBy, err);
  dim[0] = IM_MAXDIM;
  ObitInfoListPut (in->info, "BLC", OBIT_long, dim, blc, err);
  ObitInfoListPut (in->info, "TRC", OBIT_long, dim, trc, err);

  /* Open input */
  ObitImageOpen (in, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Blank fill buffer */
  ObitFArrayFill (in->image, fblank);

  /* How many input planes? */
  nplane = in->myDesc->inaxes[2];

  /* Loop over output planes */
  for (iplane=0; iplane<nplane; iplane++) {
    /* Write plane */
    ObitImageWrite (in, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* Close input */
  ObitImageClose (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* release buffer */
  in->image = ObitFArrayUnref(in->image);

} /* end ObitImageUtilBlankFill */
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
 * Interpolate selected rows from one image onto another, 
 * possibly including weighting across the image
 * Precalculated input pixel in output image may be passed
 * If doWeight, outData will be filled with the pixel values
 * interpolated from inData multiplied by a weight based on a 
 * circle defined by radius from the center; this is 1.0 
 * in the center and tapers with distance^2 to 0.0 outside
 * and the primary beam gains are written into wtData.
 * Magic value blanking supported.
 * Callable as thread
 * \param arg Pointer to InterpFuncArg argument with elements:
 * \li inDesc   Image Descriptor for input image
 * \li inData   ObitFArray with input plane pixel data
 * \li outDesc  Image Descriptor for output image
 * \li outData  ObitFArray for output plane pixel data
 * \li XPix     ObitFArray for X pixels, if NULL calculate
 * \li YPix     ObitFArray for y pixels
 * \li doWeight gboolean if TRUE, also do primary beam weighting
 * \li wtData   ObitFArray for output weight plane pixel data
 *              only used if doWeight
 * \li radius   Radius in pixels of weighting circle
 *              only used if doWeight
 * \li nZern    If>0 apply cernike corrections to position
 *              nZern is the number of terms in ZCoef
 * \li zCoef    Zernike correction coefficients
 *              only used if zCoef>0
 * \li first    First (1-rel) row in image to process this thread
 * \li last     Highest (1-rel) row in image to process this thread
 * \li ithread  thread number, <0 -> no threading
 * \li err      ObitErr Obit error stack object
 * \li thread   thread Object
 * \li Interp   ObitFInterpolate Input Image Interpolator
 * \return NULL
 */
static gpointer ThreadImageInterp (gpointer args)
{
  /* Get arguments from structure */
  InterpFuncArg *largs = (InterpFuncArg*)args;
  ObitImageDesc *inDesc = largs->inDesc;
  /* ObitFArray *inData    = largs->inData;*/
  ObitImageDesc *outDesc= largs->outDesc;
  ObitFArray *outData   = largs->outData;
  gboolean   doWeight   = largs->doWeight;
  ObitFArray *wtData    = largs->wtData;
  ObitFArray *XPix      = largs->XPixData;
  ObitFArray *YPix      = largs->YPixData;
  olong      radius     = largs->radius;
  olong      nZern      = largs->nZern;
  ofloat     *ZCoef     = largs->ZCoef;
  olong      loRow      = largs->first;
  olong      hiRow      = largs->last;
  ObitErr    *err       = largs->err;
  ObitThread *thread    = largs->thread;
  ObitFInterpolate *interp = largs->Interp;
  /* local */
  olong ix, iy, pos[2];
  ollong indx;
  ofloat *out, *outWt=NULL, rad2=0.0, dist2, irad2=0.0;
  ofloat crpix[2], wt, val, *xp=NULL, *yp=NULL, fblank =  ObitMagicF();
  ofloat inPixel[2], outPixel[2], offPixel[2];
  gboolean OK=TRUE, doPixel, sameGrid;
  gchar *routine = "ThreadImageInterp";

  /* Previous error? */
  if (err->error) goto finish;

  /* Are the grids the same? */
  sameGrid = (nZern==0) && ObitImageDescAligned(inDesc, outDesc, err);
  if (sameGrid) {
    inPixel[0] = 0.0; inPixel[1] = 0.0;
    OK = ObitImageDescCvtPixel (outDesc, inDesc, inPixel, offPixel, err);
    OK = TRUE;  /* Reset */
    /* further sanity check, must be integer */
    if (offPixel[0]>0.0) ix = (olong)(offPixel[0]+0.5);
    else                 ix = (olong)(offPixel[0]-0.5);
    sameGrid = sameGrid && (fabs(offPixel[0]-ix)<0.001);
    if (offPixel[1]>0.0) iy = (olong)(offPixel[1]+0.5);
    else                 iy = (olong)(offPixel[1]-0.5);
    sameGrid = sameGrid && (fabs(offPixel[1]-iy)<0.001);
  }
  if (err->error) {  /* Oh bother! */
    ObitThreadLock(thread);  /* Lock against other threads */
    Obit_log_error(err, OBIT_Error,"%s: Error determining alignment",
		   routine);
    ObitThreadUnlock(thread); 
    goto finish;
  }

  /* Get output aray pointer */
  pos[0] = pos[1] = 0;
  out   = ObitFArrayIndex (outData, pos);
  doPixel = XPix!=NULL;
  if (doPixel) {
    xp    = ObitFArrayIndex (XPix, pos);
    yp    = ObitFArrayIndex (YPix, pos);
  }

  /* Coordinate reference pixel of input */
  crpix[0] = inDesc->crpix[0] - inDesc->xPxOff;
  crpix[1] = inDesc->crpix[1] - inDesc->yPxOff;

  /* if weighting */
  if (doWeight) {
    /* Working version of radius */
    rad2 = radius * radius;
    irad2 = 1.0 / rad2;
    outWt = ObitFArrayIndex (wtData, pos);
  }

  /* Loop over image interpolating */
  for (iy = loRow; iy<=hiRow; iy++) { /* loop in y */
    outPixel[1] = (ofloat)iy;
    for (ix = 1; ix<=outDesc->inaxes[0]; ix++) {/* loop in x */
      outPixel[0] = (ofloat)ix;

      /* array index in out for this pixel */
      indx = (iy-1) * outDesc->inaxes[0] + (ix-1);
      
      /* Get pixel in input image*/
      if (doPixel) {   /* Precalculated? */
	inPixel[0] = xp[indx];
	inPixel[1] = yp[indx];
	OK = (inPixel[0]>=0.0) && (inPixel[1]>=0.0);
      } else {  /* Calculate  - Zernike correction? */
	if (nZern>0) { /* yes */
	  OK = ObitImageDescCvtZern (outDesc, inDesc, nZern, ZCoef, 
				     outPixel, inPixel, err);
	} else {       /* not precalculated */
	  if (sameGrid) {  /* Simple shift */
	    inPixel[0] = outPixel[0] + offPixel[0];
	    inPixel[1] = outPixel[1] + offPixel[1];
	    /* Check that in image */
	    OK = (inPixel[0]>=0.0) && (inPixel[1]>=0.0) &&
	      (inPixel[0]<inDesc->inaxes[0]) && (inPixel[0]<inDesc->inaxes[1]);
	  } else { /* Need to convert */
	  OK = ObitImageDescCvtPixel (outDesc, inDesc, outPixel, inPixel, err);
	  }
	}
	if (err->error) {
	  ObitThreadLock(thread);  /* Lock against other threads */
	  Obit_log_error(err, OBIT_Error,"%s: Error projecting pixel",
			 routine);
	  ObitThreadUnlock(thread); 
	  goto finish;
	}
      } /* End get Pixel */

      if (doWeight) { /* weighting? */
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
      } else wt = 1.0;
      
      /* interpolate */
      if (wt != fblank ) {
	val = ObitFInterpolatePixel (interp, inPixel, err);
	if (doWeight & (val != fblank )) val *= wt;
	out[indx] = val;
	if (err->error) {
	  ObitThreadLock(thread);  /* Lock against other threads */
	  Obit_log_error(err, OBIT_Error,"%s: Error interpolating pixel in %s",
			 routine, interp->name);
	  ObitThreadUnlock(thread); 
	  goto finish;
	}
	if (doWeight) outWt[indx] = wt;
      } else {
	out[indx]   = fblank;
	if (doWeight) outWt[indx] = fblank;
      }
      
    } /* end loop over x */
  } /* end loop over y */
  
  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* ThreadImageInterp */

/**
 * Threaded get input pixels in InImage for outImage.
 * Callable as thread
 * \param arg Pointer to InterpFuncArg argument with elements:
 * \li inDesc   Image Descriptor for input image
 * \li inData   ObitFArray with input plane pixel data
 * \li outDesc  Image Descriptor for output image
 * \li outData  ObitFArray for output plane pixel data
 * \li XPix     ObitFArray for X pixels
 * \li YPix     ObitFArray for y pixels
 * \li nZern    If>0 apply cernike corrections to position
 *              nZern is the number of terms in ZCoef
 * \li zCoef    Zernike correction coefficients
 *              only used if zCoef>0
 * \li first    First (1-rel) row in image to process this thread
 * \li last     Highest (1-rel) row in image to process this thread
 * \li ithread  thread number, <0 -> no threading
 * \li err      ObitErr Obit error stack object
 * \li thread   thread Object
 * \li Interp   ObitFInterpolate Input Image Interpolator
 * \return NULL
 */
static gpointer ThreadGetXYPixels (gpointer args)
{
  /* Get arguments from structure */
  InterpFuncArg *largs = (InterpFuncArg*)args;
  ObitImageDesc *inDesc = largs->inDesc;
  /* ObitFArray *inData    = largs->inData;*/
  ObitImageDesc *outDesc= largs->outDesc;
  /*ObitFArray *outData   = largs->outData;*/
  ObitFArray *XPix      = largs->XPixData;
  ObitFArray *YPix      = largs->YPixData;
  olong      nZern      = largs->nZern;
  ofloat     *ZCoef     = largs->ZCoef;
  olong      loRow      = largs->first;
  olong      hiRow      = largs->last;
  ObitErr    *err       = largs->err;
  ObitThread *thread    = largs->thread;

  /* local */
  olong ix, iy, pos[2];
  ollong indx;
  ofloat inPixel[2], outPixel[2], *xp, *yp;
  gboolean OK;
  gchar *routine = "ThreadGetXYPixels";

  /* Get output aray pointers */
  pos[0] = pos[1] = 0;
  xp    = ObitFArrayIndex (XPix, pos);
  yp    = ObitFArrayIndex (YPix, pos);

  /* Loop over image determining pixel numbers */
  for (iy = loRow; iy<=hiRow; iy++) { /* loop in y */
    outPixel[1] = (ofloat)iy;
    for (ix = 1; ix<=outDesc->inaxes[0]; ix++) {/* loop in x */
      outPixel[0] = (ofloat)ix;

     /* Get pixel in input image - Zernike correction?*/
      if (nZern>0) { /* yes */
	OK = ObitImageDescCvtZern (outDesc, inDesc, nZern, ZCoef, 
				 outPixel, inPixel, err);
      } else {       /* no */
	OK = ObitImageDescCvtPixel (outDesc, inDesc, outPixel, inPixel, err);
      }
      if (err->error) {
	ObitThreadLock(thread);  /* Lock against other threads */
	Obit_log_error(err, OBIT_Error,"%s: Error projecting pixel",
		       routine);
	ObitThreadUnlock(thread); 
	goto finish;
      }
      /* Save - array index in outout for this pixel */
      indx = (iy-1) * outDesc->inaxes[0] + (ix-1);
      if (OK) {
	xp[indx] = inPixel[0];
	yp[indx] = inPixel[1];
      } else {  /* Not in image */
	xp[indx] = -1;
	yp[indx] = -1;
      }

    } /* end loop over x */
  } /* end loop over y */
  
  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* ThreadGetXYPixels */

/**
 * Make arguments for Threaded ThreadImageInterp
 * \param thread     ObitThread object to be used for interpolator
 * \param radius     if > 0 then the radius for weighting
 * \param nZern      Number of Zernike coefficients in ZCoef
 *                   <=0 -> no Zernike corrections
 * \param ZCoef      Zernike coefficients
 * \param inDesc     input image descriptor
 * \param outDesc    output image descriptor
 * \param Interp     interpolator for input image
 *                   Cloned for multiple threads
 * \param err        Obit error stack object.
 * \param ThreadArgs[out] Created array of InterpFuncArg, 
 *                   delete with KillInterpFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
static olong MakeInterpFuncArgs (ObitThread *thread, olong radius,
				 olong nZern, ofloat *ZCoef,
				 ObitImageDesc *inDesc, ObitImageDesc *outDesc, 
				 ObitFInterpolate *Interp,
				 ObitErr *err, InterpFuncArg ***ThreadArgs)
{
  olong i, j, nThreads;

  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(thread));

  /* Initialize threadArg array */
  *ThreadArgs = g_malloc0(nThreads*sizeof(InterpFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*ThreadArgs)[i] = g_malloc0(sizeof(InterpFuncArg)); 
  for (i=0; i<nThreads; i++) {
    (*ThreadArgs)[i]->inDesc   = ObitImageDescRef(inDesc);
    (*ThreadArgs)[i]->inData   = NULL;
    (*ThreadArgs)[i]->outDesc  = ObitImageDescRef(outDesc);
    (*ThreadArgs)[i]->outData  = NULL;
    (*ThreadArgs)[i]->XPixData = NULL;
    (*ThreadArgs)[i]->YPixData = NULL;
    (*ThreadArgs)[i]->doWeight = radius>0;
    (*ThreadArgs)[i]->wtData   = NULL;
    (*ThreadArgs)[i]->radius   = radius;
    (*ThreadArgs)[i]->nZern    = nZern;
    if (nZern>0) {
      (*ThreadArgs)[i]->ZCoef = g_malloc0(nZern*sizeof(ofloat));
      for (j=0; j<nZern; j++) (*ThreadArgs)[i]->ZCoef[j] = ZCoef[j];
    } else (*ThreadArgs)[i]->ZCoef = NULL;
    (*ThreadArgs)[i]->first    = 1;
    (*ThreadArgs)[i]->last     = inDesc->inaxes[1];
    if (Interp!=NULL) {
      if (i==0) (*ThreadArgs)[i]->Interp = ObitFInterpolateRef(Interp);
      else (*ThreadArgs)[i]->Interp      = ObitFInterpolateClone(Interp, NULL);
    }
    (*ThreadArgs)[i]->ithread  = i;
    (*ThreadArgs)[i]->thread   = thread;
    (*ThreadArgs)[i]->err      = err;
  }

  return nThreads;
} /*  end MakeInterpImageArgs */

/**
 * Delete arguments for ThreadImageInterp
 * \param nargs      number of elements in args.
 * \param ThreadArgs Array of InterpFuncArg
 */
static void KillInterpFuncArgs (olong nargs, InterpFuncArg **ThreadArgs)
{
  olong i;

  if (ThreadArgs==NULL) return;
  for (i=0; i<nargs; i++) {
    if (ThreadArgs[i]) {
      if (ThreadArgs[i]->inDesc)   ObitImageDescUnref(ThreadArgs[i]->inDesc);
      if (ThreadArgs[i]->inData)   ObitFArrayUnref(ThreadArgs[i]->inData);
      if (ThreadArgs[i]->outDesc)  ObitImageDescUnref(ThreadArgs[i]->outDesc);
      if (ThreadArgs[i]->outData)  ObitFArrayUnref(ThreadArgs[i]->outData);
      if (ThreadArgs[i]->XPixData) ObitFArrayUnref(ThreadArgs[i]->XPixData);
      if (ThreadArgs[i]->YPixData) ObitFArrayUnref(ThreadArgs[i]->YPixData);
      if (ThreadArgs[i]->wtData)   ObitFArrayUnref(ThreadArgs[i]->wtData);
      if (ThreadArgs[i]->Interp)   ObitFInterpolateUnref(ThreadArgs[i]->Interp);
      if (ThreadArgs[i]->ZCoef)    g_free(ThreadArgs[i]->ZCoef);
      g_free(ThreadArgs[i]);
    }
  }
  g_free(ThreadArgs);
} /*  end KillInterpFuncArgs */

/**
 * Set Lagrangian interpolation kernal taking into account ends of the grid.
 * \param  Target Which is the desired fractional plane (1-rel)?
 * \param  npln   Number of planes
 * \param  hwidth Half width of convolution kernal
 * \param  Start  [out] first plane (1-rel) to include
 * \param  Kernal [out] convolving kernal.
 */
static void SetConvKernal3 (ofloat Target, olong npln, olong hwidth, 
			    olong *Start, ofloat *Kernal)
{
  ofloat prod, sum, xx;
  ofloat denom[10];
  olong ipos, i, j, cen, iwid;

  /* Init Lagrangian denominators for hwidth */
  iwid = 1 + (2*hwidth);
  for (j= 1; j<=iwid; j++) {
    prod = 1.0;
    for (i= 1; i<=iwid; i++) {
      if (i != j) prod = prod * (j - i);
    } 
    denom[j-1] = 1.0 / prod;
  } 

  /* fractional pixel */
  ipos = Target + 0.5;
  iwid = hwidth*2 + 1;

  /* set first pixel */
  cen = ipos - hwidth;
  cen = MAX (1, MIN (cen, (npln-iwid+1)));
  *Start = cen; /* returned version */
  /* make 0 rel */
  cen = cen - 1;

  /* set "x" at first pixel to 1.0 */
  xx = Target - cen;

  /* compute interpolating kernal */
  sum = 0.0;
  for (j= 0; j<iwid; j++) {
    prod = denom[j];
    for (i= 0; i<iwid; i++) {
      if (i != j) prod = prod * (xx - (i+1));
    } /* end i loop  */
    Kernal[j] = prod;
    sum += prod;
  } /* end j loop */

  /* Normalize to sum of 1.0 if needed */
  if (fabs(sum-1.0)>1.0e-4) {
     prod = 1.0/sum;
     for (i=0; i<iwid; i++) Kernal[i] *= prod;
  }
} /* end SetConvKernal3 */

/**
 * Fill Primary beam array
 * Callable as thread
 * \param arg Pointer to InterpFuncArg argument with elements:
 * \li thread   ObitThread object to be used for interpolator
 * \li ithread  thread number, <0 -> no threading
 * \li bs       Beam Shape object to use
 * \li outDesc  output image descriptor
 * \li outArr   Array to fill
 * \li doInvert Invert beam factor?
 * \li err      Obit error stack object.
 * \li first    First (1-rel) row in image to process this thread
 * \li last     Highest (1-rel) row in image to process this thread
 * \li thread   thread Object
 * \return NULL
 */
static gpointer ThreadImagePBCor (gpointer args)
{
  /* Get arguments from structure */
  PBCorFuncArg *largs = (PBCorFuncArg*)args;
  ObitFArray *outArr    = largs->outArr;
  ObitImageDesc *outDesc= largs->outDesc;
  ObitBeamShape *bs     = largs->bs;
  gboolean   doInvert   = largs->doInvert;
  olong      loRow      = largs->first-1;
  olong      hiRow      = largs->last-1;
  ObitErr    *err       = largs->err;
  ObitThread *thread    = largs->thread;
  /* local */
  olong ix, iy, pos[2];
  ollong indx;
  ofloat *out, fblank =  ObitMagicF();
  odouble ra, dec, dist;
  ofloat pbf, inPixel[2];
  gboolean bad=FALSE;
  gchar *routine = "ThreadImagePBCor";

  /* Previous error? */
  if (err->error) goto finish;


 /* Loop over image  */
  for (iy = loRow; iy<=hiRow; iy++) { /* loop in y */
    inPixel[1] = (ofloat)iy;
    /* Get output aray pointer */
    pos[0] = 0; pos[1] = iy;
    out = ObitFArrayIndex (outArr, pos);
    for (ix = 1; ix<=outArr->naxis[0]; ix++) {/* loop in x */
      inPixel[0] = (ofloat)ix;
      /* array index in in and out for this pixel */
      indx = (ix-1);

      /* Convert pixel to position */
      bad = 
	ObitSkyGeomWorldPos(inPixel[0], inPixel[1],
			    outDesc->crval[outDesc->jlocr], outDesc->crval[outDesc->jlocd],
			    outDesc->crpix[outDesc->jlocr], outDesc->crpix[outDesc->jlocd],
			    outDesc->cdelt[outDesc->jlocr], outDesc->cdelt[outDesc->jlocd],
			    outDesc->crota[outDesc->jlocd], &outDesc->ctype[outDesc->jlocr][4],
			    &ra, &dec);
      if (bad!=0) {  /* Oh bother! */
	ObitThreadLock(thread);  /* Lock against other threads */
	Obit_log_error(err, OBIT_Error, 
		       "%s: Error %d determining location of pixel", routine, bad);
	ObitThreadUnlock(thread); 
	goto finish;
      }
      /* Separation from pointing center */
      dist = ObitBeamShapeAngle(bs, ra, dec, 0.0);

      /* primary beam correction */
      pbf = ObitBeamShapeGainSym (bs, dist);
      if (doInvert && (fabs(pbf)>1.0e-3) && (pbf!=fblank)) pbf = 1.0/pbf;  /* Invert? */

      /* Save beam image */
      out[indx] = pbf;

    } /* end loop over x */
  } /* end loop over y */
  

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* ThreadImagePBCor */

/**
 * Make arguments for Threaded ThreadImagePBCor
 * \param thread     ObitThread object to be used for interpolator
 * \param bs         Beam Shape object to use
 * \param outDesc    output image descriptor
 * \param outArr     Array to fill
 * \param doInvert   Invert beam factor?
 * \param outDesc    output image descriptor
 * \param err        Obit error stack object.
 * \param ThreadArgs[out] Created array of PBCorFuncArg, 
 *                   delete with KillPBCorFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
static olong MakePBCorFuncArgs (ObitThread *thread, 
				ObitBeamShape *bs, ObitImageDesc *outDesc, 
				ObitFArray *outArr, gboolean doInvert,
				ObitErr *err, PBCorFuncArg ***ThreadArgs)
{
  olong i, nThreads;

  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(thread));

  /* Initialize threadArg array */
  *ThreadArgs = g_malloc0(nThreads*sizeof(PBCorFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*ThreadArgs)[i] = g_malloc0(sizeof(PBCorFuncArg)); 
  for (i=0; i<nThreads; i++) {
    (*ThreadArgs)[i]->outDesc  = ObitImageDescRef(outDesc);
    (*ThreadArgs)[i]->outArr   = ObitFArrayRef(outArr);
    (*ThreadArgs)[i]->bs       = ObitBeamShapeCopy(bs, NULL, err);
    (*ThreadArgs)[i]->ithread  = i;
    (*ThreadArgs)[i]->thread   = thread;
    (*ThreadArgs)[i]->err      = err;
    (*ThreadArgs)[i]->first    = 1;
    (*ThreadArgs)[i]->last     = outDesc->inaxes[1];
    (*ThreadArgs)[i]->doInvert = doInvert;
  }

  return nThreads;
} /*  end MakePBCorImageArgs */

/**
 * Delete arguments for ThreadImagePBCor
 * \param nargs      number of elements in args.
 * \param ThreadArgs Array of PBCorFuncArg
 */
static void KillPBCorFuncArgs (olong nargs, PBCorFuncArg **ThreadArgs)
{
  olong i;

  if (ThreadArgs==NULL) return;
  for (i=0; i<nargs; i++) {
    if (ThreadArgs[i]) {
      if (ThreadArgs[i]->outDesc)  ObitImageDescUnref(ThreadArgs[i]->outDesc);
      if (ThreadArgs[i]->outArr)   ObitFArrayUnref(ThreadArgs[i]->outArr);
      if (ThreadArgs[i]->bs)       ObitBeamShapeUnref(ThreadArgs[i]->bs);
      g_free(ThreadArgs[i]);
    }
  }
  g_free(ThreadArgs);
} /*  end KillPBCorFuncArgs */

