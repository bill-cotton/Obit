/* $Id: ObitUVCalCalibrate.h,v 1.2 2005/10/06 20:22:56 bcotton Exp $  */
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
#ifndef OBITUVCALCALIBRATE_H 
#define OBITUVCALCALIBRATE_H 

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
 * \file ObitUVCalCalibrate.h
 * ObitUVCal utilities for applying amp/phase/delay/rate calibration to
 * uv data 
 *
 * This is implemented as utility routines and a separate Structure definition
 * to avoid circular definitions. 
 * The ObitUVCal must be visible here but the structure needed for this calibration 
 * is a member of the ObitUVCal.  If it were implemented as an Obit class this would
 * lead to a circular definition which c cannot deal with.
 */

/*---------------Public functions---------------------------*/

/** Public: Init amp/phase/delay/rate Calibration */
void 
ObitUVCalCalibrateInit (ObitUVCal *in, ObitUVSel *sel, ObitUVDesc *desc, 
		     ObitErr *err);

/** Public: Calibrate data for amp/phase/delay/rate */
void ObitUVCalCalibrate (ObitUVCal *in, float time, olong ant1, 
			 olong ant2, ofloat *RP, ofloat *visIn, ObitErr *err);

/** Public: Shutdown Calibration */
void  ObitUVCalCalibrateShutdown (ObitUVCal *in, ObitErr *err);

/** Public:  Destroy structure for amp/phase/delay/rate Calibration. */
ObitUVCalCalibrateS* ObitUVCalCalibrateSUnref (ObitUVCalCalibrateS *in);

#endif /* OBITUVCALCALIBRATE_H */ 
