/* $Id:  $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009                                               */
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
/*  Define the basic components of the ObitSkyModelVMBeam structure   */
/*  This class represents sky models and their Fourier transform      */
/*  This is intended to be included in a class structure definition.  */
/**
 * \file ObitSkyModelVMBeamDef.h
 * ObitSkyModel structure members for this and any derived classes.
 */
#include "ObitSkyModelVMDef.h"  /* Parent class definitions */
/** Threshold flux density for doing high accuracy DFT model */
ofloat Threshold;
/** Current maximum residual  */
ofloat maxResid;
/** Beginning time (d) for validity of model model in VMComps */
ofloat begVMModelTime;
/** Dimension of Rgain...  */
olong dimGain;
/** Array of time/spatially variable R component gain */
ofloat *Rgain;
/** Array of time/spatially variable L component gain */
ofloat *Lgain;
/** Array of time/spatially variable Stokes Q component gain */
ofloat *Qgain;
/** Array of time/spatially variable Stokes U component gain */
ofloat *Ugain;
/** Array booleans, TRUE if corresponding antenna (0-rel) is EVLA */
gboolean *isEVLA;
/** Number of subarrays */
olong numAntList;
/** Antenna List for parallactic angle */
ObitAntennaList **AntList;
/** Current source */
ObitSource *curSource;
/** I Beam image */
ObitFullBeam *IBeam;
/** V Beam image */
ObitFullBeam *VBeam;
/** Q Beam image if doCrossPol */
ObitFullBeam *QBeam;
/** U Beam image if doCrossPol */
ObitFullBeam *UBeam;
/** Array booleans, TRUE if corresponding antenna (0-rel) is EVLA */
gboolean doCrossPol;
/** Number of Beam planes */
olong numPlane;
/** Number of UV channels */
olong numUVChann;
/** Plane for each frequency in corresponding UV data  dim numUVChann */
olong *FreqPlane;
/** Save Stokes request */
gchar saveStokes[8];
/** Save doCalSelect */
gboolean saveDoCalSelect;
/** Save doCalib */
olong saveDoCalib;
