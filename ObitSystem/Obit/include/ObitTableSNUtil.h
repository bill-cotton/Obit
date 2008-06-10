/* $Id: ObitTableSNUtil.h,v 1.3 2007/01/31 15:25:24 bcotton Exp $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2008                                          */
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
#ifndef OBITTABLESNUTIL_H 
#define OBITTABLESNUTIL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitUV.h"
#include "ObitData.h"
#include "ObitTableSN.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableSNUtil.h
 * ObitTableSN class utility routine definition.
 */

/*---------------Public functions---------------------------*/
/** Public: Copy SN tables with selection */
ObitIOCode ObitTableSNSelect (ObitUV *inUV, ObitUV *outUV, ObitErr *err);

/** Public: Invert gains in a Solution table (to revert calibration). */
ObitTableSN* ObitTableSNUtilInvert (ObitTableSN *inSN, ObitData *outData, olong *outVer, 
				    ObitErr *err);

#endif /* OBITTABLESNUTIL_H */ 
