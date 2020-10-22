/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009-2020                                          */
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
/** Number of antenna types */
olong numAntType;
/** Dimension of Rgain...  */
olong dimGain;
/** Arrays of time/spatially variable R/X component gain, real, imag */
ofloat **Rgain, **Rgaini;
/** Arrays of time/spatially variable L/Y component gain, real, imag */
ofloat **Lgain, **Lgaini;
/** Arrays of time/spatially variable Stokes RL/XY component gain, real, imag */
ofloat **RLgain, **RLgaini;
/** Arrays of time/spatially variable Stokes LR/YX component gain, real, imag */
ofloat **LRgain, **LRgaini;
/** Number of subarrays */
olong numAntList;
/** Antenna List for parallactic angle */
ObitAntennaList **AntList;
/** Antenna type index per antenna */
olong *AntType;
/** Antenna diameter (m) per type */
ofloat *Diams;
/** Current source */
ObitSource *curSource;
/**  Beam Shape */
ObitBeamShape *BeamShape;
/** R/X Beam image interpolator array */
ObitImageInterp **RXBeam;
/** L/Y Beam image interpolator array */
ObitImageInterp **LYBeam;
/** RL Beam image interpolator array  if doCrossPol */
ObitImageInterp **RLBeam;
/** LR Beam image interpolator  array if doCrossPol */
ObitImageInterp **LRBeam;
/** R/X Beam phase image interpolator array - NULL if not given */
ObitImageInterp **RXBeamPh;
/** L/Y Beam phase image interpolator array - NULL if not given **/
ObitImageInterp **LYBeamPh;
/** RL Beam phase image interpolator array if doCrossPol - NULL if not given **/
ObitImageInterp **RLBeamPh;
/** LR Beam phase image interpolator array if doCrossPol - NULL if not given **/
ObitImageInterp **LRBeamPh;
/** cross polarized corrections? */
gboolean doCrossPol;
/** Phase corrections? */
gboolean doPhase;
/** Number of Beam planes per antenna type */
olong *numPlane;
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
/** Is this model calculation in a CLEAN? */
gboolean doBeamCorClean;
