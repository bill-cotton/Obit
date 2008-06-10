/* $Id: ObitOTFGetSoln.h,v 1.7 2008/02/27 15:47:10 bcotton Exp $  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITOTFGETSOLN_H 
#define OBITOTFGETSOLN_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitOTF.h"
#include "ObitTableOTFSoln.h"
#include "ObitTableOTFCal.h"
#include "ObitImage.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFGetSoln.h
 *
 * Routines determine calibration for an ObitOTF and writes ObitTableOTFSoln
 *  tables.
 */

/*---------------Public functions---------------------------*/
/** Public: Determine offset calibration from a residual  data. */
ObitTableOTFSoln* ObitOTFGetSolnCal (ObitOTF *in, ObitOTF *out, ObitErr *err);

/** Public: Determine Gain calibration from a residual  data. */
ObitTableOTFSoln* ObitOTFGetSolnGain (ObitOTF *in, ObitOTF *out, ObitErr *err);

/** Public: Offset calibration from filtered residual  data. */
ObitTableOTFSoln* ObitOTFGetSolnFilter (ObitOTF *in, ObitOTF *out, ObitErr *err);

/** Public: Fit polynomial baseline to median filtered residual OTF data. */
ObitTableOTFSoln* ObitOTFGetSolnPolyBL (ObitOTF *in, ObitOTF *out, ObitErr *err);

/** Public: Determine multibeam baseline calibration from  data. */
ObitTableOTFSoln* ObitOTFGetSolnMBBase (ObitOTF *in, ObitOTF *out, ObitErr *err);

/** Public: Determine median value of the cal from a set of data. */
ofloat ObitOTFGetSolnAvgCal (gint n, gboolean *isCal, ofloat *data);

/** Public: Determine instrumental calibration. */
ObitTableOTFSoln* ObitOTFGetInstCal (ObitOTF *in, ObitOTF *out, ObitErr *err);

/** Public: Determine Penn Array type gain calibration from data. */
ObitTableOTFSoln* ObitOTFGetSolnPARGain (ObitOTF *in, ObitOTF *out, ObitErr *err);

/** Public: Fill Soln table with pointing corrections. */
ObitTableOTFSoln* ObitOTFGetSolnPointTab (ObitOTF *inOTF, ObitOTF *outOTF, 
					  ObitErr *err);

/** Public: Create Dummy Cal table. */
ObitTableOTFCal* ObitOTFGetDummyCal (ObitOTF *in, ObitOTF *out, olong ver, 
				     oint ncoef, ObitErr *err);

/** Public: Flag data on basis of comparison of statistics of model, residuals. */
void ObitOTFGetSolnFlag (ObitOTF *inOTF, ObitImage *model, 
			 ObitOTF *outOTF, olong FGVer, ObitErr *err);
#endif /* OBITOTFGETSOLN_H */ 
