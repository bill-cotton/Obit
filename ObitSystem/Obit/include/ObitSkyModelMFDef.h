/* $Id: ObitSkyModelMFDef.h 119 2009-07-17 16:02:35Z bill.cotton $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010                                               */
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
/*  Define the basic components of the ObitSkyModel structure         */
/*  This class represents sky models and their Fourier transform      */
/*  This is intended to be included in a class structure definition.   */
/**
 * \file ObitSkyModelDef.h
 * ObitSkyModelMF structure members for this and any derived classes.
 */
#include "ObitSkyModelDef.h"  /* Parent class definitions */
/** Number of spectral elements  */
olong nSpec;
/** Array of component model with tabulated spectra as rows in an FArray */
ObitFArray *specComps;
/** Array of nSpec central frequencies in Hz */
odouble *specFreq;
/** Array of CC spectrum index (0-rel) per input vis channel */
olong *specIndex;
/** Array of image planes of single image */
ObitFArray **planes;
/** Array of Fourier transforms of planes */
ObitCArray **FTplanes;
/** Array of Interpolators for UV grids */
ObitCInterpolate **myInterps;
