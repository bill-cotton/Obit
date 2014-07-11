/* $Id$ */
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
/*  Define the basic components of the ObitSkyModelVMPoln structure   */
/*  This class represents sky models and their Fourier transform      */
/*  This is intended to be included in a class structure definition.  */
/**
 * \file ObitSkyModelVMPolnDef.h
 * ObitSkyModel structure members for this and any derived classes.
 */
#include "ObitSkyModelVMDef.h"  /* Parent class definitions */
/** Number of Stokes to use of (I,Q,U,V) - typically 3 or 4 */
olong nModelStokes;
/** Number of coarse spectral planes in Stokes I,Q,U,V models */
olong nSpecI, nSpecQ, nSpecU, nSpecV;
/** Prior spectral index for Stokes I, Q, U, V */
ofloat AlphaI, AlphaQ, AlphaU, AlphaV;
/** Prior spectral index reference frequency (Hz) for Stokes I, Q, U, V */
odouble AlphaIRefF, AlphaQRefF, AlphaURefF, AlphaVRefF;
/** Stokes I, Q, U, V image mosaic */
ObitImageMosaic *mosaicI, *mosaicQ, *mosaicU, *mosaicV;
/** Array of Stokes I, Q, U, V components as rows in an FArray */
ObitFArray *VMIComps, *VMQComps, *VMUComps, *VMVComps;
/** Array of MF frequency bins (Hz), Stokes I, Q, U, V */
odouble *specFreqI, *specFreqQ, *specFreqU, *specFreqV;
/** Index for MF Frequency bin per channel/IF, Stokes I, Q, U, V */
olong *specIndexI, *specIndexQ, *specIndexU, *specIndexV;
/** Number of actual components in VMIComps... */
olong numIComp, numQComp, numUComp, numVComp;
/** List of beginning I, Q, U, V  component per image in mosaic (1-rel) */
olong *startIComp, *startQComp, *startUComp, *startVComp;
/** List of highest I, Q, U, V  component per image in mosaic (1-rel) */
olong *endIComp, *endQComp, *endUComp, *endVComp;
/** Array booleans, TRUE if corresponding antenna (0-rel) is EVLA */
gboolean *isEVLA;
/** Number of antennas */
olong numAnt;
/** Number of subarrays */
olong numAntList;
/** Antenna List for parallactic angle */
ObitAntennaList **AntList;
/** Current source */
ObitSource *curSource;
/** Number of UV channels */
olong numUVChan;
/** Plane for each frequency in corresponding UV data  dim numUVChann */
olong *FreqPlane;
/** Save Stokes request */
gchar saveStokes[8];
/** Save doCalSelect */
gboolean saveDoCalSelect;
/** Save doCalib */
olong saveDoCalib;
/** PD Table */
olong PDVer;
/** PD Reference Antenna */
olong PDrefAnt;
/** Beginning time (d) for validity of parallactic angles */
ofloat begVMModelTime;
/** Is Circular feed array */
gboolean isCirc;
/** Instrumental polarization type */
ObitUVPolCalType polType;
/** Array of frequency, antenna RS, RSc inst. poln term  arrays Circular feeds
 [numUVChann][numAnt] */
dcomplex **RS, **RSc;
/** Array of frequency, antenna RD, RDc arrays Circular feeds */
dcomplex **RD, **RDc;
/** Array of frequency, antenna LS, LSc arrays Circular feeds */
dcomplex **LS, **LSc;
/** Array of frequency, antenna LD, RLc arrays Circular feeds */
dcomplex **LD, **LDc;
/** Array of frequency cross pol ref ant factors - Circular feeds */
dcomplex *PPLR, *PPRL;
/** Array of frequency, antenna CX, CXc inst. poln term  arrays Linear feeds
 [numUVChann][numAnt] */
dcomplex **CX, **CXc;
/** Array of frequency, antenna SX, SXc arrays  Linear feeds */
dcomplex **SX, **SXc;
/** Array of frequency, antenna CY, CYc arrays  Linear feeds */
dcomplex **CY, **CYc;
/** Array of frequency, antenna SY, SYc arrays Linearr feeds */
dcomplex **SY, **SYc;
/** Array of frequency phasedifference factors factors - Linear feeds */
dcomplex *PPXY, *PPYX;
