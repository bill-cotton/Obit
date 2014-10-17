/* $Id:  $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012                                               */
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
#ifndef OBITTABLEFSUTIL_H 
#define OBITTABLEFSUTIL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitTableFS.h"
#include "ObitTableVL.h"
#include "ObitImage.h"
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableFSUtil.h
 * ObitTableFS class utility routine definition.
 */

/*---------------Public functions---------------------------*/
/** Public: Append one FS table to another */
void ObitTableFSAppend (ObitTableFS *in, ObitTableFS *out, ObitErr *err);

/** Public: Index a FS table */
void ObitTableFSIndex (ObitTableFS *in, ObitErr *err);

/** Public: Merge overlapping components */
void ObitTableFSMerge (ObitTableFS *in, ObitErr *err);

/** Public: Select significant components */
void ObitTableFSSelect (ObitTableFS *in, ObitTableFS *out, ObitErr *err);

/** Public: Remove entries from a given field */
void ObitTableFSPurge (ObitTableFS *in, gchar *field, ObitErr *err);

/** Public: Remove redundant entries */
void ObitTableFSRedun (ObitTableFS *in, ObitTableFS *out, ObitErr *err);

/** Public: Apply final calibration and error analysis */
void ObitTableFSCal (ObitTableFS *in, ObitErr *err);

/** Public: Write human readable version of an FS table to a FILE */
void ObitTableFSPrint (ObitTableFS *in, ObitImage *image, FILE  *prtFile, 
		       ObitErr *err);

/** Public: Loop through image extracting spectrum */
void  ObitTableFSGetSpectrum(ObitTableFS *in, ObitImage *im, ObitErr *err);

/** Public: Convert VL to FS */
ObitTableFS* ObitTableFSUtilVL2FS (ObitTableVL *in, ObitData *data, olong FSver, 
				   ObitErr *err);

/** Public: Filter and fit spectra */
void ObitTableFSFiltVel(ObitTableFS *inFS, ObitImage *im, ObitTableFS *outFS, 
			ObitErr *err);

#endif /* OBITTABLEFSUTIL_H */ 
