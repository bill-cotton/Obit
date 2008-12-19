/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                          */
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
/*  Define the basic components of the ObitSkyModel ClassInfo structure */
/* This is intended to be included in a classInfo structure definition  */
#include "ObitClassDef.h"  /* Parent class ClassInfo definition file */
/** Function pointer to Create SkyModel object from description in an ObitInfoList. */
ObitSkyModelFromInfoFP ObitSkyModelFromInfo;
/** Function pointer to Constructor. */
ObitSkyModelCreateFP ObitSkyModelCreate;
/** Function pointer to Initializer. */
ObitSkyModelInitModFP ObitSkyModelInitMod;
/** Function pointer to Finalizer. */
ObitSkyModelShutDownModFP ObitSkyModelShutDownMod;
/** Function pointer to initialize model for pass in time through data  */
ObitSkyModelInitModelFP ObitSkyModelInitModel;
/** Function pointer to  Load specified image and plane */
ObitSkyModelLoadFP ObitSkyModelLoad;
/** Function pointer to Subtract model from an ObitUV */
ObitSkyModelSubUVFP ObitSkyModelSubUV;
/** Function pointer to  Divide model into an ObitUV */
ObitSkyModelDivUVFP ObitSkyModelDivUV;
/** Function pointer to Calculate Fourier transform of model for 
    current buffer in an ObitUV*/
ObitSkyModelFTFP ObitSkyModelFT;
/** Function pointer to Load point model */
ObitSkyModelLoadPointFP ObitSkyModelLoadPoint;
/** Function pointer to Load Components model */
ObitSkyModelLoadCompsFP ObitSkyModelLoadComps;
/** Function pointer to Grid Components model */
ObitSkyModelGridCompsFP ObitSkyModelGridComps;
/** Function pointer to Load image model */
ObitSkyModelLoadImageFP ObitSkyModelLoadImage;
/** Function pointer to DFT FT */
ObitSkyModelFTDFTFP ObitSkyModelFTDFT;
/** Function pointer to Grid FT */
ObitSkyModelFTGridFP ObitSkyModelFTGrid;
/** Function pointer to  Sum flux in Clean Model */
ObitSkyModelSumFP ObitSkyModelSum;
/** Function pointer to Compress CC Tables */
ObitSkyModelCompressCCFP ObitSkyModelCompressCC;
/** Function pointer to Get input parameters */
ObitSkyModelGetInputFP ObitSkyModelGetInput;
/** Function pointer to Decide model method */
ObitSkyModelChoseFP ObitSkyModelChose;
/** Function pointer to Fill in data selection values */
ObitSkyModelSetSelectFP ObitSkyModelSetSelect;
/** Function pointer to Decide next block of channels for PB corr */
ObitSkyModelsetPBChansFP ObitSkyModelsetPBChans;
/** Function pointer to ObitTableCC with possible PB corrections  */
ObitSkyModelgetPBCCTabFP ObitSkyModelgetPBCCTab;
/** Function pointer to fill in->plane with image and possibly PB corrected */
ObitSkyModelgetPBImageFP ObitSkyModelgetPBImage;
/** Function pointer to  Grid/FT components */
ObitSkyModelGridFTCompsFP ObitSkyModelGridFTComps;
/** Function pointer to Load Grid component */
ObitSkyModelLoadGridCompsFP ObitSkyModelLoadGridComps;
/** Function pointer to FT image array  */
ObitSkyModelFTImageFP ObitSkyModelFTImage;
/** Function pointer to AddField  */
ObitSkyModelAddFieldFP ObitSkyModelAddField;
/** Function pointer to Extract information about underlying structures to ObitInfoList */
ObitSkyModelGetInfoFP ObitSkyModelGetInfo;
