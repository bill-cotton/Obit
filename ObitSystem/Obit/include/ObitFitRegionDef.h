/* $Id: ObitFitRegionDef.h,v 1.1 2006/05/12 14:05:38 bcotton Exp $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006                                               */
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
/*  Define the basic components of the ObitFitRegion structure        */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitFitRegionDef.h
 * ObitFitRegion structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/** bottom left corner in selected region of image (0-rel) */
olong corner[2];
/** dimension of region */
olong dim[2];
/** peak in region */
ofloat peak;
/** peak in region residual after model subtraction */
ofloat peakResid;
/** RMS residual */
ofloat RMSResid;
/** Sum of pixel values in residual */
ofloat fluxResid;
/** Number of models */
olong nmodel;
/** Array of Models */
ObitFitModel **models;
