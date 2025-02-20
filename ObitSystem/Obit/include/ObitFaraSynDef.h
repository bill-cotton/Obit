/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2025                                               */
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
/*  Define the basic components of the ObitFaraSyn structure          */
/*  This is intended to be included in a class structure definition.  */
/**
 * \file ObitFaraSynDef.h
 * ObitFaraSyn structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class instance definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** Number of Number of Lambda^2 planes */
olong nLamb2;
/** Array of wavelength^2 (m^2) in inQFArray, inUFArray */
odouble *lamb2;
/** Reference wavelength^2 for synthesis/analysis */
odouble refLamb2;
/** Array of QPol FArrays in Lambda^2 */
ObitFArray **inQFArrays;
/** Array of UPol FArrays in Lambda^2 */
ObitFArray **inUFArrays;
/** Array of output arrays */
ObitFArray **outFArrays;
/** Weights for each Q/U plane */
ofloat *plnWt;
/**  RMS of each Q plane */
ofloat *QRMS;
/**  RMS of each U plane */
ofloat *URMS;
/** Median and sigma QRMS */
ofloat qMedian, qSigma;
/** Median and sigma URMS */
ofloat uMedian, uSigma;
/** Number of Faraday depth planes */
olong nOut;
/** Min. Faraday depth (rad/m^2) */
ofloat minRMSyn;
/** Max. Faraday depth (rad/m^2) */
ofloat maxRMSyn;
/** Increment in Faraday depth (rad/m^2) */
ofloat delRMSyn;
/** CArray Of Faraday depth */
ObitCArray **RMPlanes;
/** Output Amp image for 3D Faraday Synthesis */
ObitImage* outAImage;
/** Output Phase image for 3D Faraday Synthesis */
ObitImage* outPImage;
/** Output Cube for 2D  Faraday Analysis */
ObitImage* outFitRM;
/** If true write Amplitudes to outAImage
    Also phases to outPImage if non NULL */
gboolean doWrite;
/** SNR for Q and U pixels - Faraday Analysis */
ofloat minQUSNR;
/** Min. fraction of planes included - Faraday Analysis */
ofloat minFrac;
/** If true do error analysis - Faraday Analysis */
gboolean doError;
/** If true do max RM synthesis - Faraday Analysis */
gboolean doRMSyn;
/** Max. Chi^2 for search - Faraday Analysis */
ofloat maxChi2;
