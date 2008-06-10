/* $Id: ObitIoN2SolNTable.h,v 1.3 2007/09/07 12:31:39 bcotton Exp $ */
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
#ifndef OBITION2SOLNTABLE_H 
#define OBITION2SOLNTABLE_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitUV.h"
#include "ObitTableNI.h"
#include "ObitTableSN.h"
#include "ObitInfoList.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIoN2SolNTable.h
 * ObitIoN2SolNTable utility definition.
 *
 * This utility package converts an IoNospheric model into a SolutioN table.
 * There are no objects of this class.
 */


/*---------------Public functions---------------------------*/
/** Public: Convert a NI table to an SN table. */
ObitTableSN* ObitIoN2SolNTableConvert (ObitUV *inUV, ObitTableNI *NITable,
				       ObitTableSN *outSN, ofloat shift[2], 
				       ObitErr *err);

#endif /* ION2SOLNTABLE_H */ 
