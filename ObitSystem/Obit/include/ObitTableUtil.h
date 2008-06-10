/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                          */
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
#ifndef OBITTABLEUTIL_H 
#define OBITTABLEUTIL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitTable.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableUtil.h
 * ObitTable class utility routine definition.
 */

/*---------------Public functions---------------------------*/
/**  Public: Sort a table by one key */
ObitIOCode ObitTableUtilSort (ObitTable *in, gchar *colName, gboolean desc,
			      ObitErr *err);
/**  Public: Sort a table by abs value of one key */
ObitIOCode ObitTableUtilAbsSort (ObitTable *in, gchar *colName, gboolean desc,
				 ObitErr *err);
/**  Public: Sort a table by otwo numeric keys */
ObitIOCode ObitTableUtilSort2f (ObitTable *in, 
				gchar *colName1, olong cell1, gboolean desc1,
				gchar *colName2, olong cell2, gboolean desc2, 
				ObitErr *err);

#endif /* OBITTABLEUTIL_H */ 
