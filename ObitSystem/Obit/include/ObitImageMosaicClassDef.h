/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2011                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
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
/*  Define the basic components of the ObitImageMosaic ClassInfo structure  */
/* This is intended to be included in a classInfo structure definition*/
#include "ObitClassDef.h"  /* Parent class ClassInfo definition file */
/** Function pointer to Constructor from ObitInfo. */
ObitImageMosaicFromInfoFP ObitImageMosaicFromInfo;
/** Function pointer to Extract information about underlying structures 
    to ObitInfoList. */
ObitImageMosaicGetInfoFP ObitImageMosaicGetInfo;
/** Function pointer to Zap specified image. */
ObitImageMosaicZapImageFP ObitImageMosaicZapImage;
/** Function pointer to Return specified image */
ObitImageMosaicGetImageFP ObitImageMosaicGetImage;
/** Function pointer to Set specified image */
ObitImageMosaicSetImageFP ObitImageMosaicSetImage;
/** Function pointer to Return RMS pixel value of  image */
ObitImageMosaicGetImageRMSFP ObitImageMosaicGetImageRMS;
/** Function pointer to ObitImageMosaicGetFullImage */
ObitImageMosaicGetFullImageFP ObitImageMosaicGetFullImage;
/** Function pointer to Set  Full Field  image */
ObitImageMosaicSetFullImageFP ObitImageMosaicSetFullImage;
/** Function pointer to Set underlying files */
ObitImageMosaicSetFilesFP ObitImageMosaicSetFiles;
/** Function pointer to Create Mosaic from uv data */
ObitImageMosaicCreateFP ObitImageMosaicCreate;
/** Function pointer to Define parameters of images */
ObitImageMosaicDefineFP ObitImageMosaicDefine;
/** Function pointer to Flatten tiles onto full field image  */
ObitImageMosaicFlattenFP ObitImageMosaicFlatten;
/** Function pointer to Give field of view */
ObitImageMosaicFOVFP ObitImageMosaicFOV;
/** Function pointer to Give max. field of view in current model */
ObitImageMosaicMaxFOVFP ObitImageMosaicMaxFOV;
/** Function pointer to Get max summed CC and determine offset from nearest pixel */
ObitImageMosaicMaxCCFP ObitImageMosaicMaxCC;
/** Function pointer to Give combined Clean Components table */
ObitImageMosaicCombineCCFP ObitImageMosaicCombineCC;
/** Function pointer to Zero selected CC entries */
ObitImageMosaicFlagCCFP ObitImageMosaicFlagCC;
/** Function pointer to Add field to mosaic */
ObitImageMosaicAddFieldFP ObitImageMosaicAddField;
/** Function pointer to Generate a mosaic for peeling */
ObitImageMosaicMaxFieldFP ObitImageMosaicMaxField;
/** Function pointer to Concatenate Image CC tables onto the FullField Image */
ObitImageMosaicCopyCCFP ObitImageMosaicCopyCC;
/* Private functions for derived classes */
/** Function pointer to  Cover specified field of view */
FlyEyeFP FlyEye;
/**  Function pointer to Add field to list */
AddFieldFP AddField;
/** Function pointer to Lookup outliers in catalog */
AddOutlierFP AddOutlier;
/** Function pointer to add beam tapers */
AddTapersFP AddTapers;
