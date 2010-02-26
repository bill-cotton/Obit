/*   $Id: ObitUVGridMFDef.h 76 2009-02-04 14:51:56Z bill.cotton $ */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitUV structure               */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitUVGridMFDef.h
 * ObitUVGridMF structure members for this and any derived classes.
 */
#include "ObitUVGridDef.h"  /* Parent class definitions */
/** Maximum order of the imaging 
    Spectral index only = 1, plus curvature = 2 */
olong maxOrder;
/** Number of coarse frequency planes */
olong nSpec;
/** Arrays of start and finish IFs (0-rel), per coarse channel */
olong *BIFSpec, *EIFSpec;
/** Arrays of start and finish Channels (0-rel), per coarse channel */
olong *BChanSpec, *EChanSpec;
/** Array of accumulation grids, per coarse channel */
ObitCArray **grids;
/** Array of Beam normalization factors, per coarse channel */
ofloat *BeamNorms;
/** Arrays gridding tapers per uv data channel for forcing the beams */
ofloat *sigma1, *sigma2, *sigma3;
