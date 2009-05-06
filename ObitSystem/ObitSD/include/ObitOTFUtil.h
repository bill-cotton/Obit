/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2009                                          */
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
#ifndef OBITOTFUTIL_H 
#define OBITOTFUTIL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitIO.h"
#include "ObitImage.h"
#include "ObitOTF.h"
#include "ObitOTFArrayGeom.h"
#include "ObitOTFSkyModel.h"
#include "ObitOTFGrid.h"
#include "ObitFInterpolate.h"
#include "ObitTableCC.h"
#include "ObitTableOTFSoln.h"

/*-------- Obit:  Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFUtil.h
 * Utility routines for the ObitOTF class.
 *
 * \section ObitOTFUtilparameters Control Parameters
 * The imaging control parameters are passed through the info object 
 * on the OTF data, these control both the output image files and the 
 * processing parameters.
 * OTF Data selection/calibration/editing control
 * \li  "doCalSelect" OBIT_bool (1,1,1) Select/calibrate/edit data?
 * \li  "doCalib" OBIT_int (1,1,1) >0 -> calibrate,
 * \li  "gainUse" OBIT_int (1,1,1) SN/CL table version number, 0-> use highest
 * \li  "flagVer" OBIT_int (1,1,1) Flag table version, 0-> use highest, <0-> none
 * \li  "Stokes" OBIT_string (4,1,1) Selected output Stokes parameters:
 *               "I", "V", " " -> "I"
 * \li  "BChan" OBIT_int (1,1,1) First spectral channel selected. [def all]
 * \li  "EChan" OBIT_int (1,1,1) Highest spectral channel selected. [def all]
 * \li  "Targets" OBIT_string (?,?,1) Target names selected. [def all]
 * \li  "timeRange" OBIT_float (2,1,1) Selected timerange in days. [def all]
 * \li  "Scans" OBIT_int (2,1,1) Lowest and highest selected scan numbers. [def all]
 * \li  "Feeds" OBIT_int (?,1,1) a list of selected feed numbers, [def all.]
 * 
 * Gridding/imaging parameters
 * \li "minWt"     OBIT_float (1,1,1) Minimum summed gridding convolution weight 
 *                 as a fraction of the maximum [def 0.01]
 */

/*---------------Public functions---------------------------*/
/** Public: Subtract an image (ObitFarray) from an ObitOTF. */
void ObitOTFUtilSubImage(ObitOTF *inOTF, ObitOTF *outOTF, ObitFArray *image, 
			 ObitImageDesc *desc, ObitErr *err);

/** Public: Replace data with an image (ObitFarray) model. */
void ObitOTFUtilModelImage(ObitOTF *inOTF, ObitOTF *outOTF, ObitFArray *image, 
			   ObitImageDesc *desc, ObitErr *err);

/** Public: Scale an OTF data set. */
void ObitOTFUtilScale(ObitOTF *inOTF, ObitOTF *outOTF, ofloat scale, ofloat offset,
		      ObitErr *err);

/** Public: Scale and add Gaussian noise to an OTF data set. */
void ObitOTFUtilNoise(ObitOTF *inOTF, ObitOTF *outOTF, ofloat scale, ofloat offset,
		      ofloat sigma, ObitErr *err);

/** Public: Replace data with a sky model in a buffer of data. */
void ObitOTFUtilModelImageBuff (ObitOTF *in, ObitFInterpolate *image, ofloat factor, 
				ObitErr *err);

/** Public: Subtract a sky model from a buffer of data. */
void ObitOTFUtilSubSkyModelBuff (ObitOTF *in, ObitOTFSkyModel *sky, ofloat factor);

/** Public: Create an Image Object from an OTF */
ObitImage* ObitOTFUtilCreateImage (ObitOTF *inOTF, ObitErr *err);

/** Public: Convert an OTF to an image Object */
void ObitOTFUtilMakeImage (ObitOTF *inOTF,  ObitImage *outImage, gboolean doBeam,
			   ObitImage *Beam, ObitImage *Wt, ObitErr *err);

/** Public:Index an OTF */
void ObitOTFUtilIndex (ObitOTF *inOTF, ObitErr *err);

/** Public: Get on-off differences in a nodding scan */
void ObitOTFUtilDiffNod (ObitOTF *inOTF, olong scan, ObitErr *err);

/** Public:Create image cube */
void 
ObitOTFUtilMakeCube (ObitImageDesc *inDesc, ObitOTFDesc *OTFDesc, 
		     ObitImageDesc *outDesc, 
		     gchar *Stokes, olong bchan, olong echan, olong incr, ObitErr *err);

/** Public:  Utility to convolve CCs with a beam */
ObitFArray* ObitOTFUtilConvBeam (ObitTableCC *CCTab, ObitImage *Beam, 
				 ObitFArray *Template, ObitErr *err);
/** Public:  Residual calibration */
ObitTableOTFSoln* ObitOTFUtilResidCal (ObitOTF *inOTF, ObitOTF *outOTF, 
				       ObitImage *model, gboolean doModel,
				       ObitImage *PSF, ObitErr *err);
#endif /* OBITOTFUTIL_H */ 
