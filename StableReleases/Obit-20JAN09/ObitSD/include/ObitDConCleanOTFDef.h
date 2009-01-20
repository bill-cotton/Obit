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
/*  Define the basic components of the ObitDConCleanOTF structure        */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitDConCleanOTFDef.h
 * ObitDConCleanOTF structure members for this and any derived classes.
 */
#include "ObitDConCleanDef.h"  /* Parent class definitions */
/** CLEAN Restoring beam size FWHM in pixels */
ofloat cleanSize;
/** Dirty Image - use instead of mosaic */
ObitImage *dirty;
/** Beam Image - Raw instrumental response */
ObitImage *beam;
/** Clean Image */
ObitImage *clean;
/** Gridding Weight Image */
ObitImage *weight;
/** Restore image when done? */
gboolean doRestore;
/** Scale CCs by ration of dirty to CLEAN beam areas? */
gboolean doScaleCC;
