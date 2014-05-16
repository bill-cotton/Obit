/* $Id$        */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITCUDASKYMODEL_H 
#define OBITCUDASKYMODEL_H 
#include "ObitCUDASkyModelInfoDef.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCUDASkyModel.h
 *
 * Define interface to (primitive) CUDA routines
 * Sky models are calculated using a GPU.
 * 
 */

/*---------------Public functions---------------------------*/
/* Public: Initialize DFT Model */
void ObitCUDASkyModelDFTInit (GPUInfo *gpuInfo, GPUVisInfo *visInfo, GPUModelInfo *modelInfo);

/* Public: Set DFT sky model Model */
void ObitCUDASkyModelDFTSetMod (GPUInfo *gpuInfo, GPUVisInfo *visInfo, GPUModelInfo *modelInfo);

/* Public: Calculate DFT Model */
void ObitCUDASkyModelDFTCalc (GPUInfo *gpuInfo, GPUVisInfo *visInfo, GPUModelInfo *modelInfo);

/* Public: Shutdown DFT Model */
void ObitCUDASkyModelDFTShutdown (GPUInfo *gpuinfo, GPUVisInfo   *visInfo, GPUModelInfo *modelInfo);
#endif /* OBITFCUDASKYMODEL_H */ 
