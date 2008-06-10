/* $Id: ObitDConCleanPxListDef.h,v 1.6 2006/07/17 16:00:05 bcotton Exp $ */
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
/*;  Correspondence concerning Obit should be addressed as follows:   */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitDConCleanPxList structure  */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitDConCleanBmHistDef.h
 * ObitDConCleanPxList structure members for this and any derived classes.
 */
#include "ObitDef.h"  /* Parent class (Obit) definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** Image mosaic being deconvolved */
ObitImageMosaic *mosaic;
/** Image Clean window definition */
ObitDConCleanWindow *window;
/** Image Beam patch data,
    The beam patch array is square and 2*patch +1 x 2*patch +1
    in size.  The center pixel is the center of the dirty beam.*/
ObitFArray *BeamPatch;
/** Minimum flux density to load */
ofloat minFluxLoad;
/** Maximum abs. residual at end of CLEAN */
ofloat maxResid;
/** min. minor cycle flux for auto Window feature */
ofloat autoWinFlux ;
/** Total flux density CLEANed */
ofloat totalFlux;
/** Number of fields */
olong nfield;
/** Plane being processed, 1-rel indices of axes 3-? */
olong plane[IM_MAXDIM-2];
/** Maximum number of iterations */
olong niter;
/** current number of iterations */
olong currentIter;
/** current number of components per field */
olong *iterField;
/** CC table version per field */
olong *CCver;
/** current summed flux density per field */
ofloat *fluxField;
/** Gaussian size (circular) per field, (deg.) */
ofloat *circGaus;
/** CLEAN gain per field */
ofloat *gain;
/** Minimum flux density per field */
ofloat *minFlux;
/** Depth factor per field */
ofloat *factor;
/** Array of CCTables for writing */
ObitTableCC **CCTable;
/** maximum number of pixels */
olong maxPixel;
/** number of pixels */
olong nPixel;
/** pixel x coordinate (0-rel) */
olong *pixelX;
/** pixel y coordinate (0-rel) */
olong *pixelY;
/** pixel field */
gshort *pixelFld;
/** pixel flux density */
ofloat *pixelFlux;
/** message level for deconvolution progress messages, 
    0=none, 1=summary, 2=normal, higher numbers for diagnostics */
olong prtLv;

