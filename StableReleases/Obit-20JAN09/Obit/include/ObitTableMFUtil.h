/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006                                               */
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
#ifndef OBITTABLEMFUTIL_H 
#define OBITTABLEMFUTIL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitTableMF.h"
#include "ObitTableVL.h"
#include "ObitFitRegionList.h"
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableMFUtil.h
 * ObitTableMF class utility routine definition.
 */

/*---------------Public functions---------------------------*/
/** Public: Convert an ObitFitRegionList to MF table entries */
void ObitTableMFRegions2MF (ObitTableMF *in, ObitFitRegionList *regList, 
			    ObitImage *image, ObitErr *err);

/** Public: Convert an MF table to a VL Table */
void ObitTableMF2VL (ObitTableMF *in, ObitTableVL *out, 
		     ObitImage *image, ObitErr *err);

/** Public: Write human readable version of an MF table to a FILE */
void ObitTableMFPrint (ObitTableMF *in, ObitImage *image, FILE  *prtFile, 
		       ObitErr *err);
#endif /* OBITTABLEMFUTIL_H */ 
