/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2008                                          */
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
#ifndef OBITTABLEVLUTIL_H 
#define OBITTABLEVLUTIL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitTableVL.h"
#include "ObitTableVZ.h"
#include "ObitFitRegionList.h"
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableVLUtil.h
 * ObitTableVL class utility routine definition.
 */

/*---------------Public functions---------------------------*/
/** Public: Append one VL table to another */
void ObitTableVLAppend (ObitTableVL *in, ObitTableVL *out, ObitErr *err);

/** Public: Index a VL table */
void ObitTableVLIndex (ObitTableVL *in, ObitErr *err);

/** Public: Merge overlapping components */
void ObitTableVLMerge (ObitTableVL *in, ObitErr *err);

/** Public: Select significant components */
void ObitTableVLSelect (ObitTableVL *in, ObitTableVL *out, ObitErr *err);

/** Public: Remove entries from a given field */
void ObitTableVLPurge (ObitTableVL *in, gchar *field, ObitErr *err);

/** Public: Remove redundant entries */
void ObitTableVLRedun (ObitTableVL *in, ObitTableVL *out, ObitErr *err);

/** Public: Apply final calibration and error analysis */
void ObitTableVLCal (ObitTableVL *in, ObitErr *err);

/** Public: Write human readable version of an VL table to a FILE */
void ObitTableVLPrint (ObitTableVL *in, ObitImage *image, FILE  *prtFile, 
		       ObitErr *err);

/** Public: Convert a VL table to a VZ (short) table */
ObitTableVZ* ObitTableVL2VZ (ObitTableVL *in, ObitData *data, ObitErr *err);

/** Public: Select entries in a VZ (short) table */
ObitTableVZ* ObitTableVZSel (ObitTableVZ *in,  ObitData *data, ObitErr *err);
#endif /* OBITTABLEVLUTIL_H */ 
