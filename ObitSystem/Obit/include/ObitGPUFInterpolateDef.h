/* $Id: ObitGPUSkyModelDef.h 490 2014-08-15 18:21:51Z bill.cotton $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014                                               */
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
/*  Define the basic components of the ObitGPUSkyModel structure      */
/*  This is intended to be included in a class structure definition.  */
/**
 * \file ObitGPUSkyModelDef.h
 * ObitGPUSkyModel structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/* Needs definitions in ObitCUDASkyModelInfoDef.h */
/** Dimension of plane in myArray */
int nx, ny;
/** Half width of interpolation kernal */
int hwidth;
/** Object in CUDA space */
CUDAFInterpolate *FInterpolate;
/** Input Array */
ObitGPUFArray *inArray;
/** X/Y pixel arrays */
ObitGPUFArray *xArray, *yArray;
/** output locked memory buffer */
ofloat *buffer;
