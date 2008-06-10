/* $Id: ObitDConCleanDef.h,v 1.4 2005/04/19 11:46:30 bcotton Exp $ */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitDConClean structure        */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitDConCleanDef.h
 * ObitDCon structure members for this and any derived classes.
 */
#include "ObitDConDef.h"  /* Parent class definitions */
/** CLEAN window list */
ObitDConCleanWindow *window;
/** Current Beam patch */
ObitFArray *BeamPatch;
/** Current Pixel list */
ObitDConCleanPxList *Pixels;
/** Current Beam histogram */
ObitDConCleanBmHist *BeamHist;
/** Current Pixel histogram */
ObitDConCleanPxHist *PixelHist;
/** Current field number 1-rel */
olong currentField;
/** Beam patch halfwidth size in pixels */
olong beamPatchSize;
/** Minimum beam patch halfwidth in pixels */
olong minPatchSize;
/** Minimum flux density to load in this CLEAN cycle*/
ofloat minFluxLoad;
/** Maximum number of iterations */
olong niter;
/** Maximum number of residuals */
olong maxPixel;
/** Number to skip in loading Pixel list */
olong numberSkip;
/** Restoring beam, values in deg*/
ofloat bmaj, bmin, bpa;
/** Number of fields in arrays */
olong nfield;
/** CLEAN gain per field */
ofloat *gain;
/** Minimum flux density per field */
ofloat *minFlux;
/** Depth factor per field */
ofloat *factor;
/** CC table version same for each field */
olong CCver;
/** Max Abs windowed residual per field, -1=> uninitialized */
ofloat *maxAbsRes;
/** Average windowed residual per field, -1=> uninitialized */
ofloat *avgRes;
/** auto Window feature requested? */
gboolean autoWindow;
/** min. minor cycle flux for auto Window feature */
ofloat autoWinFlux ;
