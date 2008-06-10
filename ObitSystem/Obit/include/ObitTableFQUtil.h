/* $Id$ */
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
#ifndef OBITTABLEFQUTIL_H 
#define OBITTABLEFQUTIL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitUV.h"
#include "ObitTableFQ.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableFQUtil.h
 * ObitTableFQ class utility routine definition.
 */

/*---------------Public functions---------------------------*/
/** Public: Read information for an FQid  */
ObitIOCode ObitTableFQGetInfo (ObitTableFQ *in, oint fqid, oint *nif,
			       odouble **freqOff, oint **sideBand,
			       ofloat **chBandw, ObitErr *err);

/** Public: Write Read information for an FQid */
ObitIOCode ObitTableFQPutInfo (ObitTableFQ *in, oint fqid, oint nif,
			       odouble *freqOff, oint *sideBand,
			       ofloat *chBandw, ObitErr *err);

/** Public: Copy FQ table with selection */
ObitIOCode ObitTableFQSelect (ObitUV *inUV, ObitUV *outUV, odouble *SouIFOff,
			      odouble SouBW, ObitErr *err);
#endif /* OBITTABLEFQUTIL_H */ 
