/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004                                               */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitCInterpolate structure     */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitCInterpolateDef.h
 * ObitInterpolate structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** Image descriptor to give relation between pixels and coordinates */
ObitImageDesc *myDesc;
/** Array to be interpolated */
ObitCArray *myArray;
/** Pointer to data in myArray */
ofloat *array;
/** Dimension of plane in myArray */
olong nx, ny;
/** Number of conjugate (neg V) Columns in myArray  */
olong numConjCol;
/** Table of interpolation kernals, every kernalSpace of a pixel */
ObitFArray *myKernal;
/** Spacing of kernal tabulation */
ofloat kernalSpace;
/** Number of tabulations in myKernal per cell = 1/kernalSpace */
olong numKTab;
/** Half width of interpolation kernal */
olong hwidth;
/** Target "X" Pixel */
ofloat xPixel;
/** Target "Y" Pixel */
ofloat yPixel;
/** Convolving "X" start pixel number in convolution */
olong xStart;
/** Convolving "Y" start pixel number in convolution */
olong yStart;
/** Number of  "X" terms in convolution */
olong xNterm;
/** Number of  "Y" terms in convolution  */
olong yNterm;
/** "X" Convolving kernal (pointer in kernalTable) */
ofloat *xKernal;
/** "Y" Convolving kernal (pointer in kernalTable)*/
ofloat *yKernal;
/** Reciprocals of Lagrangian denominators */
ofloat denom[40];
