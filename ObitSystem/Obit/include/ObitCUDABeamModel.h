/* $Id$    */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2026                                               */
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
#ifndef OBITCUDABEAMMODEL_H 
#define OBITCUDABEAMMODEL_H
#include "ObitCUDABeamModelInfoDef.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCUDABeamModel.h
 *
 * CUDA implemention of skyModel with Beam corrections
 * Define interface to (primitive) CUDA routines
 * Sky models are calculated using a GPU.
 * Real GPU/CUDA versions are called if HAVE_GPU==1, else stubbed versions
 * 
 */

/*---------------Public functions---------------------------*/
#if HAVE_GPU==1  /* GPU? Real versions */
/* Public: Initialize DFT Model */
#ifdef __cplusplus
extern "C" {
#endif
void ObitCUDABeamModelDFTInit (GPUBeamInfo *gpuInfo);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
/* Public: Set DFT sky model Model */
void ObitCUDABeamModelDFTSetMod (GPUBeamInfo *gpuInfo);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
/* Public: Calculate DFT Model */
  void ObitCUDABeamModelDFTCalc (GPUBeamInfo *gpuInfo, int doJones, int doUnit,
			       float parAng, long visRange[2]);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
/* Public: Free device CUDAFInterpolate */
void CUDAFInterpolateFree (CUDAFInterpolate* in);
#ifdef __cplusplus
}
#endif

/* Public: Shutdown DFT Model */
#ifdef __cplusplus
extern "C" {
#endif
void ObitCUDABeamModelDFTShutdown (GPUBeamInfo *gpuinfo);
#ifdef __cplusplus
}
#endif

#else  /* No GPU - stubb */
/* Public: Initialize DFT Model */
void ObitCUDABeamModelDFTInit (GPUBeamInfo *gpuInfo)
{
} /* end ObitCUDABeamModelDFTInit */

/* Public: Set DFT sky model Model */
void ObitCUDABeamModelDFTSetMod (GPUBeamInfo *gpuInfo)
{
} /* end  ObitCUDABeamModelDFTSetMod  */

/* Public: Calculate DFT Model */
void ObitCUDABeamModelDFTCalc (GPUBeamInfo *gpuInfo,
			      int doJones, float parAng,
			      long visRange[2])
{
} /* end  ObitCUDABeamModelDFTCalc */

/* Public: Shutdown DFT Model */
void ObitCUDABeamModelDFTShutdown (GPUBeamInfo *gpuinfo)
{
} /* end ObitCUDABeamModelDFTShutdown */
#endif /* HAVE_GPU */
#endif /* OBITFCUDABEAMMODEL_H */ 
