/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2021                                               */
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
/*  Define the basic components of the ObitGPUGrid ClassInfo structure */
/* This is intended to be included in a classInfo structure definition  */
#include "ObitClassDef.h"  /* Parent class ClassInfo definition file */
/** Function pointer to Constructor. */
ObitGPUGridCreateFP ObitGPUGridCreate;
/** Function pointer to copy gridding data to GPU (CUDA) structures. */
ObitGPUGridSetGPUStructFP ObitGPUGridSetGPUStruct;
/** Function pointer to setuo gridding structures for GPU*/
ObitGPUGridSetGPUFP ObitGPUGridSetGPU;
/** Function pointer to initialize GPU*/
ObitGPUGridInitGPUFP ObitGPUGridInitGPU;
/** Function pointer to copy UV data to GPU */
ObitGPUGrid2GPUFP ObitGPUGrid2GPU;
/** Function pointer to Grid UV data */
ObitGPUGridGridFP ObitGPUGridGrid;
/** Function pointer to Flip UV data */
ObitGPUGridFlipFP ObitGPUGridFlip;
/** Function pointer to  copy gridded data to host */
ObitGPUGrid2HostFP ObitGPUGrid2Host;
/** Function pointer to shutdown */
ObitGPUGridShutdownFP ObitGPUGridShutdown;
