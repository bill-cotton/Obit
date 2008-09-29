/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                          */
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
 * ObitSkyModel structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** I/O status */
ObitIOStatus myStatus;
/** Image array */
ObitImageMosaic *mosaic;
/** Image plane of single image */
ObitFArray *plane;
/** Fourier transform of plane */
ObitCArray *FTplane;
/** Number of conjugate (neg U) columns in  FTplane */
olong numConjCol;
/** Array of component model as rows in an FArray */
ObitFArray *comps;
/** Interpolator for UV grid */
ObitCInterpolate *myInterp;
/** Model type */
ObitSkyModelType modelType;
/** Component model type */
ObitSkyModelCompType modType;
/** Model calculation mode for components */
ObitSkyModelMode modelMode;
/** Current model calculation mode */
ObitSkyModelMode currentMode;
/** List of AIPSCC table versions per image in mosaic 
 there are mosaic->numberImages of these */
olong *CCver;
/** List of beginning component per image in mosaic (1-rel) */
olong *startComp;
/** List of highest component per image in mosaic (1-rel) */
olong *endComp;
/** Factor to multiply times model */
ofloat factor;
/** Minimum flux density model or pixel */
ofloat minFlux;
/** Factor to multiply times second Stokes of model */
ofloat stokFactor;
/** Point model flux density (Jy) */
ofloat pointFlux;
/** Point, x (ra), y (dec) offset in deg. */
ofloat pointXOff, pointYOff;
/** Other (non-point)model components:
    major_axis (deg),  minor_axis (deg),  position_angle (deg),
    type (ObitSkyModelCompType as gint), spectral terms;
 */
ofloat pointParms[10];
/** Antennna diameter for rel. PB corrections*/
ofloat antSize;
/** Apply 3D imaging corrections */
gboolean do3D;
/** Divide model into data? */
gboolean doDivide;
/** Replace data with model? */
gboolean doReplace;
/** Make relative Primary Beam corrections? */
gboolean doPBCor;
/** Selected start channel[1-rel] and number */
olong startChannel, numberChannel;
/** Selected start IF [1-rel] and number  */
olong startIF, numberIF;
/** Selected start rel. PB correction channel[1-rel] and number */
olong startChannelPB, numberChannelPB;
/** Selected start rel. PB correction IF [1-rel] and number  */
olong startIFPB, numberIFPB;
/** Number of frequency channels for PB correction */
olong nfreqPB;
/** Reference frequency for this block of channels for PB corrections */
odouble PBFreq;
/** Selected Stokes */
gchar stokes[5];
/** Selected start Poln [1-rel] and number  */
olong startPoln, numberPoln;
/** True if need to multiply the FT by sqrt(-1) before applying */
gboolean doFlip;
/** True if only positive flux components are to be used */
gboolean noNeg;
/** Minimum absolute component flux to use in DFT */
ofloat minDFT;
/** Maximum absolute component flux to use in Gridded model */
ofloat maxGrid;
/** Something to do for DFT model */
gboolean doDFT;
/** Something to do for Grid model  */
gboolean doGrid;
/** message level for deconvolution progress messages, 
    0=none, 1=summary, 2=normal, higher numbers for diagnostics */
olong prtLv;
/** Number of spectral terms  */
olong nSpecTerm;
/** Number of threads (elements in threadArgs)  */
olong nThreads;
/** Array of FT Function structures  */
gpointer **threadArgs;
/** DFT Fourier transform routine  */
ObitThreadFunc DFTFunc;
/** Gridded Fourier transform routine  */
ObitThreadFunc GridFunc;
