/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2008                                          */
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
/*  Define the basic components of the ObitDConCleanVis structure     */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitDConCleanVisDef.h
 * ObitDConCleanVis structure members for this and any derived classes.
 */
#include "ObitDConCleanDef.h"  /* Parent class definitions */
/** UVImager to create images from UV data */
ObitUVImager *imager;
/** Sky Model for subtracting Clean model from data */
ObitSkyModel *skyModel;
/** Model calculation mode for components */
ObitSkyModelMode modelMode;
/** Restore image when done? */
gboolean doRestore;
/** Cross Restore images when done? */
gboolean doXRestore;
/** Flatten image when done? */
gboolean doFlatten;
/** Weight data before Cleaning? */
gboolean doWeight;
/** Do beams need to be remade before cleaning? */
gboolean doBeam;
/** Consider recentering components? */
gboolean doRecenter;
/** Peak in an image window encountered */
ofloat peakFlux;
/** Quality (desirability for next CLEAN) measure per field */
ofloat *quality;
/** CLEAN component level to reuse, <0 none  */
ofloat reuseFlux;
/** autoCenter min flux density   */
ofloat autoCen;
/** Display server */
ObitDisplay *display;
/** Copy of calibrated/weighted data if doing SDI clean */
ObitUV *SDIdata;
