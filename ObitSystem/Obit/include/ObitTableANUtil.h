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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITTABLEANUTIL_H 
#define OBITTABLEANUTIL_H 

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"
#include "ObitUV.h"
#include "ObitAntennaList.h"
#include "ObitTableAN.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableANUtil.h
 * ObitTableAN class utility routine definition.
 */

/*---------------Public functions---------------------------*/
/** Public: Get Highest antenna number, reference frequency, date */
ObitIOCode ObitTableANGetInfo (ObitTableAN *in, oint *numAnt, odouble *refFreq, 
			       gchar* refDate, ObitErr *err);

/** Public: Read an antenna list into a ObitAntennaList */
ObitAntennaList* ObitTableANGetList (ObitTableAN *in, ObitErr *err);

/** Public: Copy AN tables with selection */
ObitIOCode ObitTableANSelect (ObitUV *inUV, ObitUV *outUV, ObitErr *err);
#endif /* OBITTABLEANUTIL_H */ 
