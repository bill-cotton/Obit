/* $Id: ObitUVCalSelect.h,v 1.3 2007/08/31 17:24:49 bcotton Exp $ */
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
#ifndef OBITUVCALSELECT_H 
#define OBITUVCALSELECT_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitUVDesc.h"
#include "ObitUVSel.h"
#include "ObitUVCal.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVCalSelect.h
 * ObitUVCal utilities for selecting data and translating Stokes.
 *
 * This is implemented as utility routines and a separate Structure definition
 * to avoid circular definitions. 
 * The ObitUVCal must be visible here but the structure needed for this calibration 
 * is a member of the ObitUVCal.  If it were implemented as an Obit class this would
 * lead to a circular definition which c cannot deal with.
 */

/*--------------Class definitions-------------------------------------*/

/** Public: Init data selection */
void 
ObitUVCalSelectInit (ObitUVCal *in, ObitUVSel *sel, ObitUVDesc *inDesc, 
		     ObitUVDesc *outDesc, ObitErr *err);

/** Public: translate/select data */
gboolean ObitUVCalSelect (ObitUVCal *in, ofloat *visIn, ofloat *RP, ofloat *visOut, ObitErr *err);

#endif /* OBITUVCALSELECT_H */ 
