/* $Id: ObitUVCalBandpass.h,v 1.2 2005/10/06 20:22:56 bcotton Exp $  */
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
#ifndef OBITUVCALBANDPASS_H 
#define OBITUVCALBANDPASS_H 

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"
#include "ObitUVDesc.h"
#include "ObitUVSel.h"
#include "ObitUVCal.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVCalBandpass.h
 * ObitUVCal utilities for applying bandpass calibration to uv data 
 *
 * This is implemented as utility routines and a separate Structure definition
 * to avoid circular definitions. 
 * The ObitUVCal must be visible here but the structure needed for this calibration 
 * is a member of the ObitUVCal.  If it were implemented as an Obit class this would
 * lead to a circular definition which c cannot deal with.
 */

/*---------------Public functions---------------------------*/

/** Public: Init bandpass calibration */
void 
ObitUVCalBandpassInit (ObitUVCal *in, ObitUVSel *sel, ObitUVDesc *desc, 
		     ObitErr *err);

/** Public: Apply bandpass calibration */
void ObitUVCalBandpass (ObitUVCal *in, float time, olong ant1, 
			 olong ant2, ofloat *RP, ofloat *visIn, ObitErr *err);

/** Public: Shutdown Calibration */
void  ObitUVCalBandpassShutdown (ObitUVCal *in, ObitErr *err);

/** Public:  Destroy bandpass calibration structure . */
ObitUVCalBandpassS* ObitUVCalBandpassSUnref (ObitUVCalBandpassS *in);

#endif /* OBITUVCALBANDPASS_H */ 
