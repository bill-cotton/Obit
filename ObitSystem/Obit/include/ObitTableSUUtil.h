/* $Id: ObitTableSUUtil.h,v 1.5 2007/08/31 17:24:48 bcotton Exp $ */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITTABLESUUTIL_H 
#define OBITTABLESUUTIL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitTableSU.h"
#include "ObitSourceList.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableSUUtil.h
 * ObitTableSU class utility routine definition.
 */

/*---------------Public functions---------------------------*/
/** Public: Lookup a list of sources */
ObitIOCode ObitTableSULookup (ObitTableSU *in, gint32 *dim, gchar *inlist, 
			      olong Qual, gchar souCode[5], olong *outlist, 
			      gboolean *select, olong *Number, ObitErr *err);

/** Public: Read a source list into a ObitSourceList */
ObitSourceList* ObitTableSUGetList (ObitTableSU *in, ObitErr *err);

#endif /* OBITTABLESUUTIL_H */ 
