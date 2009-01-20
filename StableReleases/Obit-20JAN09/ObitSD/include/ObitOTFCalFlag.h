/* $Id$                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003                                               */
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
#ifndef OBITOTFCALFLAG_H 
#define OBITOTFCALFLAG_H 

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"
#include "ObitOTFDesc.h"
#include "ObitOTFSel.h"
#include "ObitOTFCal.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFCalFlag.h
 * ObitOTFCal utilities for applying flagging to uv data 
 *
 * This is implemented as utility routines and a separate Structure definition
 * to avoid circular definitions. 
 * The ObitOTFCal must be visible here but the structure needed for this calibration 
 * is a member of the ObitOTFCal.  If it were implemented as an Obit class this would
 * lead to a circular definition which c cannot deal with.
 */

/*---------------Public functions---------------------------*/

/** Public: Init Flagging */
void 
ObitOTFCalFlagInit (ObitOTFCal *in, ObitOTFSel *sel, ObitOTFDesc *desc, 
		     ObitErr *err);

/** Public: Apply Flagging  */
void ObitOTFCalFlag (ObitOTFCal *in, float time, ofloat *recIn, ObitErr *err);

/** Public: Shutdown Flagging  */
void  ObitOTFCalFlagShutdown (ObitOTFCal *in, ObitErr *err);

/** Public:  Destroy Flagging structure. */
ObitOTFCalFlagS* ObitOTFCalFlagSUnref (ObitOTFCalFlagS *in);

#endif /* OBITOTFCALFLAG_H */ 
