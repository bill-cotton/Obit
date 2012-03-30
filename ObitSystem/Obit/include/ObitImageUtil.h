/* $Id$       */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2012                                          */
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
#ifndef OBITIMAGEUTIL_H 
#define OBITIMAGEUTIL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitInfoList.h"
#include "ObitUV.h"
#include "ObitImage.h"
#include "ObitImageMF.h"
#include "ObitImageWB.h"
#include "ObitImageDesc.h"
#include "ObitUVDesc.h"
#include "ObitFArray.h"
#include "ObitUVGrid.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageUtil.h
 * ObitImageUtil module definition.
 *
 * This utility module contains functions for manipulating images and 
 * imaging uv data.
 */

/*----------------- Macroes ---------------------------*/

/** 
 * Convenience Macro to set imaging parameters on an #ObitUV.
 * Sets values on ObitInfoList on input object.
 * May fail for literal constants.
 *\li in   = ObitUV to specify parameters for.
 *\li nVisPIO Number of visibilities per read/write
 *\li nChAvg Number of spectral channels to average per image 
 *           plane. (OBIT_int)
 *\li rotate Image rotation angle (same for all).(OBIT_float)
 *           0=> image all. Ignored if IF axis larger than 1.
 *\li NField Number of fields, the following all are arrays 
 *           of this dimension.
 *\li nx     Number of "X" (RA) pixels in image. (OBIT_int)
 *\li nxBeam Number of "X" (RA) pixels in beam.  (OBIT_int)
 *\li ny     Number of "Y" (Dec) pixels in image.(OBIT_int)
 *\li nyBeam Number of "Y" (Dec) pixels in beam. (OBIT_int)
 *\li xCells Cell spacing in X in asec. (OBIT_float)
 *\li yCells Cell spacing in Y in asec. (OBIT_float)
 *\li xShift Shift in X (after rotation) in deg. (OBIT_float)
 *\li yShift Shift in Y (after rotation) in deg. (OBIT_float)
 *\li err = ObitErr to receive error messages.
 */
#define ObitImageUtilSet(in,nVisPIO,nChAvg,rotate,nfield,nx,nxBeam,ny,nyBeam,xCells,yCells,xShift,yShift,err)  G_STMT_START{ \
       in->info->dim[0]=1; in->info->dim[1]=1; in->info->dim[2]=1;  \
       in->info->dim[3]=1; in->info->dim[4]=1;                      \
       in->info->work[0]= nChAvg;  in->info->work[1]= nVisPIO;      \
       in->info->fwork[0]= rotate;                                  \
       ObitInfoListPut (in->info, "nVisPIO", OBIT_long,              \
		  in->info->dim, (gpointer)&in->info->work[1], err);\
       ObitInfoListPut (in->info, "nChAvg", OBIT_long,               \
		  in->info->dim, (gpointer)&in->info->work[0], err);\
       ObitInfoListPut (in->info, "rotate", OBIT_float,             \
		  in->info->dim, (gpointer)&in->info->fwork[0],err);\
       in->info->dim[0] = nfield;                                   \
       ObitInfoListPut (in->info, "nx", OBIT_long, in->info->dim,    \
		 (gpointer)&nx, err);                               \
       ObitInfoListPut (in->info, "nxBeam", OBIT_long, in->info->dim,\
		 (gpointer)&nxBeam, err);                           \
       ObitInfoListPut (in->info, "ny", OBIT_long, in->info->dim,    \
		 (gpointer)&ny, err);                               \
       ObitInfoListPut (in->info, "nyBeam", OBIT_long, in->info->dim,\
		 (gpointer)&nyBeam, err);                           \
       ObitInfoListPut (in->info, "xCells",OBIT_float,in->info->dim,\
		 (gpointer)xCells, err);                            \
       ObitInfoListPut (in->info, "yCells",OBIT_float,in->info->dim,\
		 (gpointer)yCells, err);                            \
       ObitInfoListPut (in->info, "xShift",OBIT_float,in->info->dim,\
		 (gpointer)xShift, err);                            \
       ObitInfoListPut (in->info, "yShift",OBIT_float,in->info->dim,\
		 (gpointer)yShift, err);                            \
     }G_STMT_END  

/*---------------Public functions---------------------------*/
/** Public: Create an ObitImage from uv data. */
ObitImage* ObitImageUtilCreateImage (ObitUV *inUV, olong fieldNo,
			       gboolean doBeam, ObitErr *err);

/** Public: Fill an image with an image made from from uv data. */
void ObitImageUtilMakeImage (ObitUV *inUV, ObitImage *outImage, 
			     olong channel, gboolean doBeam, 
			     gboolean doWeight, ObitErr *err);

/** Public: Parallel fill images with those image made from from uv data. */
void ObitImageUtilMakeImagePar (ObitUV *inUV, olong nPar, ObitImage **outImage, 
			     gboolean doBeam, gboolean doWeight, ObitErr *err);

/** Public: Make an image from from uv data, info in ObitInfoList. */
void ObitImageUtilMakeImageFileInfo (ObitInfoList *inList,ObitErr *err);

/** Public: Convert an ObitUVDesc to an ObitImageDesc */
void 
ObitImageUtilUV2ImageDesc(ObitUVDesc *UVDesc, ObitImageDesc *imageDesc,
			  gboolean doGrid, olong nchavg);

/** Public: Interpolate pixels in one image to another. */
void ObitImageUtilInterpolateImage (ObitImage *inImage, ObitImage *outImage, 
				    olong *inPlane, olong *outPlane,
				    olong hwidth, ObitErr *err);

/** Public: Interpolate pixels in one image to another with Zernike corrections */
void 
ObitImageUtilInterpolateImageZern (ObitImage *inImage, ObitImage *outImage, 
				   olong *inPlane, olong *outPlane,
				   olong hwidth, olong nZern, ofloat *ZCoef, 
				   ObitErr *err);

/** Public: Interpolate pixels in one image to another giving output weighting. */
void ObitImageUtilInterpolateWeight (ObitImage *inImage, ObitImage *outImage, 
				     ObitImage *outWeight, gboolean memOnly,
				     olong radius, olong *inPlane, olong *outPlane,
				     olong hwidth, ObitErr *err);

/** Public: Correct (divide) an image by the primary beam pattern of another. */
void ObitImageUtilPBCorr (ObitImage *inImage, ObitImage *pntImage, ObitImage *outImage, 
			   olong *inPlane, olong *outPlane, ofloat antSize, ObitErr *err);

/** Public: Multiply an image by the primary beam pattern of another. */
void ObitImageUtilPBApply (ObitImage *inImage, ObitImage *pntImage, ObitImage *outImage, 
			   olong *inPlane, olong *outPlane, ofloat antSize, ObitErr *err);

/** Public: Fill image with the primary beam pattern */
void ObitImageUtilPBImage (ObitImage *pntImage, ObitImage *outImage, 
			   olong *outPlane, ofloat antSize, ofloat minGain, 
			   ObitErr *err);

/** Public: determine imaging parameters from UVW Extrema */
void ObitImageUtilImagParm (ofloat MaxBL, ofloat MaxW,
			    ofloat *Cells, ofloat *Radius);

/** Public: Create a FITS image from an ObitFArray */
ObitImage* ObitImageUtilArray2Image (gchar *fileName, olong disk, 
				     ObitFArray *inArray, ObitErr *err);

/** Public: Quantize an ObitImage and write to a FITS image */
ObitImage* ObitImageUtilQuanFITS (ObitImage *inImage, gchar *fileName, 
				  olong disk, ObitErr *err);

/** Public: Define an Image freq cube descriptor from a single plane and uv Descriptors */
void ObitImageUtilMakeCube (ObitImageDesc *inDesc, ObitUVDesc *uvDesc, 
			    ObitImageDesc *outDesc, 
			    gchar *Stokes, olong bchan, olong echan, olong incr, ObitErr *err);

/** Public: Insert a plane from an image into a cube */
void ObitImageUtilInsertPlane (ObitImage *in, ObitImage *out, olong *plane, 
			       ObitErr *err);

/** Public: Insert multiple planes from image in starting at plane in out. */
void ObitImageUtilInsertCube (ObitImage *in, ObitImage *out, olong *plane, 
			      olong axExp, ObitErr *err);

/** Public: Flux weighted velocity image from Cube. */
void ObitImageUtilVel (ObitImage *inImage, ObitImage *outImage, ObitErr *err);

/** Public: Copy with selection by pixel increment. */
void ObitImageUtilSelCopy (ObitImage *inImage, ObitImage *outImage, ObitErr *err);

/** Public: Filter out of band noise. */
void ObitImageUtilUVFilter (ObitImage *inImage, ObitImage *outImage, ofloat radius, 
			    ObitErr *err);

/** Public: Write ObitFArray as a FITS image. */
ObitImage* ObitImageUtilFArray2FITS (ObitFArray *array, 
				     gchar *FITSFile, olong FITSdisk,
				     ObitImageDesc *desc, ObitErr *err);

/** Public: Interpolate an image to a grid with a fixed shift */
void ObitImageUtilShift (ObitImage *inImage,  ObitImage *outImage, 
			 ofloat *shift, ObitErr *err);

/** Public: Interpolate an MF image to a grid with a fixed shift */
void ObitImageUtilMFShift (ObitImage *inImage,  ObitImage *outImage, 
			 ofloat *shift, ObitErr *err);

/** Public: Give parallel imaging buffer size */
olong ObitImageUtilBufSize (ObitUV *inU);

/** Public: Set Image header 2D shift parameters */
void ObitImageUtilTwoDShift (ObitUVDesc *UVDesc, ObitImageDesc *imageDesc,
			     gboolean onGrid);

/** Convert image with TSpec to Spec model type */
void ObitImageUtilT2Spec  (ObitImage *inImage, ObitImage **outImage, 
			   olong nTerm, olong *inCCVer, olong *outCCVer,
			   olong startComp, olong endComp, ObitErr *err);
/** Fit beam size to dirty beam */
void ObitImageUtilFitBeam (ObitImage *beam, ObitErr *err);

#endif /* OBITIMAGEUTIL_H */ 
