/* $Id$                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003                                               */
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
/*  Define the basic components of the ObitFInterpolate structure     */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitFInterpolateDef.h
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
ObitFArray *myArray;
/** Pointer to data in myArray */
ofloat *array;
/** Dimension of plane in myArray */
olong nx, ny;
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
/** "X" Convolving kernal */
ofloat xKernal[10];
/** "Y" Convolving kernal */
ofloat yKernal[10];
/** Reciprocals of Lagrangian denominators */
ofloat denom[10];
