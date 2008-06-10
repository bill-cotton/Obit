/* $Id$ */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitSkyModelIon structure      */
/*  This class represents sky models and their Fourier transform      */
/*  This is intended to be included in a class structure definition.  */
/**
 * \file ObitSkyModelDef.h
 * ObitSkyModel structure members for this and any derived classes.
 */
#include "ObitSkyModelDef.h"  /* Parent class definitions */
/** Array of time/spatially variable component model as rows in an FArray */
ObitFArray *VMComps;
/** Current time (d) for model in VMComps */
ofloat curVMModelTime;
/** End time (d) for validity of model model in VMComps */
ofloat endVMModelTime;
/** Model update interval (d) */
ofloat updateInterval;
/** Number of actual components in model */
olong numComp;
