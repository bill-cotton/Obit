/* $Id: ObitDConCleanPxListWBDef.h 128 2009-09-23 14:48:29Z bill.cotton $ */
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
/*;  Correspondence concerning Obit should be addressed as follows:   */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitDConCleanPxListWB structure*/
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitDConCleanPxListWBDef.h
 * ObitDConCleanPxListWB structure members for this and any derived classes.
 */
#include "ObitDConCleanPxListDef.h"  /* Parent class (Obit) definitions */
/* Max. order number of spectral beams */
olong order;
/** Current order of the imaging 
    Spectral index only = 1, plus curvature = 2 */
olong curOrder;
/** Array of Image Beam00 patch data (dirty beam convolved with dirty beam),
    Per field, the beam patch array is square and 2*patch+1 x 2*patch +1
    in size.  The center pixel is the center of the beam.    */
ObitFArray **BeamPatch00;
/**  Dirty beam convolved with spectral beam 1   */
ObitFArray **BeamPatch01;
/**  Dirty beam convolved with spectral beam 2   */
ObitFArray **BeamPatch02;
/** spectral beam 1 convolved with  Dirty beam  */
ObitFArray **BeamPatch10;
/** spectral beam 1 convolved with  spectral beam 1 */
ObitFArray **BeamPatch11;
/** spectral beam 1 convolved with spectral beam 2  */
ObitFArray **BeamPatch12;
/** spectral beam 2 convolved with  Dirty beam  */
ObitFArray **BeamPatch20;
/** spectral beam 2 convolved with  spectral beam 1 */
ObitFArray **BeamPatch21;
/** spectral beam 2 convolved with spectral beam 2  */
ObitFArray **BeamPatch22;
/** pixel flux density convolved with spectral beam 1 */
ofloat *pixelFlux1;
/** pixel flux density convolved with spectral beam 2 */
ofloat *pixelFlux2;
/** pixel flux density convolved with spectral beam 3 */
ofloat *pixelFlux3;
/** Reference frequency */
odouble refFreq;
/** Frequency at low end of the band */
odouble loFreq;
/** Frequency at high end of the band */
odouble hiFreq;
