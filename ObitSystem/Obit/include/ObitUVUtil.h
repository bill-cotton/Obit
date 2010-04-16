/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2009                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;                                                                   */
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
#ifndef OBITUVUTIL_H 
#define OBITUVUTIL_H 

#include "ObitErr.h"
#include "ObitUV.h"
#include "ObitSourceList.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVUtil.h
 * ObitUVUtil module definition.
 *
 * This utility module contains functions which operate on ObitUV
 * data sets.
 */

/*---------------Public functions---------------------------*/
/** Public: Get UVW Extrema */
void ObitUVUtilUVWExtrema (ObitUV *inUV, ofloat *MaxBL, ofloat *MaxW,
			   ObitErr *err);

/** Public: Make zeroed copy of an ObitUV */
ObitUV* ObitUVUtilCopyZero (ObitUV *inUV, gboolean scratch, ObitUV *outUV,
			    ObitErr *err);

/** Public: Divide the visibilities in one ObitUVData by another */
void ObitUVUtilVisDivide (ObitUV *inUV1, ObitUV *inUV2, ObitUV *outUV, 
			  ObitErr *err);

/** Public: Subtract the visibilities in one ObitUVData from another 
 * \relates ObitUV */
void ObitUVUtilVisSub (ObitUV *inUV1, ObitUV *inUV2, ObitUV *outUV, 
		       ObitErr *err);

/** Public: Compare the visibilities in one ObitUVData with another */
ofloat ObitUVUtilVisCompare (ObitUV *inUV1, ObitUV *inUV2, ObitErr *err);

/** Public:  Index a uv data */
void ObitUVUtilIndex (ObitUV *inUV, ObitErr *err);

/** Public: Get list of selected sources  */
ObitSourceList* ObitUVUtilWhichSources (ObitUV *inUV, ObitErr *err);

/** Public: Average a data set in frequency */
ObitUV* ObitUVUtilAvgF (ObitUV *inUV, gboolean scratch, ObitUV *outUV,
			ObitErr *err);

/** Public: Average a data set in time */
ObitUV* ObitUVUtilAvgT (ObitUV *inUV, gboolean scratch, ObitUV *outUV,
			ObitErr *err);

/** Public: Average a data set in time and/or frequency */
ObitUV* ObitUVUtilBlAvgTF (ObitUV *inUV, gboolean scratch, ObitUV *outUV,
			   ObitErr *err);

/** Public: Count good data by time segment  */
ObitInfoList* ObitUVUtilCount (ObitUV *inUV, ofloat timeInt, ObitErr *err);

/** Public: Copy data by channel to multiple output UV  */
void ObitUVUtilSplitCh (ObitUV *inUV, olong nOut, ObitUV **outUV, 
			ObitErr *err);

/** Public: Add Gaussian noise to a UV data set */
void ObitUVUtilNoise(ObitUV *inUV, ObitUV *outUV, ofloat scale, ofloat sigma, 
		     ObitErr *err);

/** Public: Add a flag entry */
ObitIOCode ObitUVUtilFlag(ObitUV *inUV, ObitErr *err);

/** Public: Calculate visibility uvw */
void ObitUVUtilUVW(const ofloat b[3], odouble dec, ofloat ha, ofloat uvw[3]);

/** Public: Append the contents of one UV data onto another */
void ObitUVUtilAppend(ObitUV *inUV, ObitUV *outUV, ObitErr *err);

/** Public: How many channels can I average? */
olong ObitUVUtilNchAvg(ObitUV *inUV, ofloat maxFact, ofloat FOV, ObitErr *err);

#endif /* OBITIUVUTIL_H */ 
