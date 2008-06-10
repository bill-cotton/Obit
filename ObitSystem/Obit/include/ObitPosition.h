/* $Id$    */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 1996,1997-2008                                     */
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
#ifndef OBITPOSITION_H 
#define OBITPOSITION_H 

#include "ObitImageDesc.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPosition.h
 * ObitPosition utilities
 *
 * Utility routines related to celestial position
 */

/*---------------Public functions---------------------------*/
/** Public: Routine to determine accurate position for pixel coordinates */
olong ObitPositionWorldPos(ofloat pixel[2], ObitImageDesc *desc,
			  odouble coord[2]);

/** Public: Routine to determine accurate pixel coordinates for an RA and Dec
 * offsets from the reference position. */
olong ObitPositionXYpix(odouble coord[2], ObitImageDesc *desc,
		       ofloat pixel[2]);

/** Public: Routine to determine accurate position for pixel coordinates from 
 * offsets from the reference position. */
olong ObitPositionWorldPosLM(odouble doff[2], ObitImageDesc *desc,
			    odouble coord[2]);

/** Public: Routine to determine accurate coordinate offsets for an RA and Dec
 * offsets from the reference position.*/
olong ObitPositionXYpixLM(odouble coord[2], ObitImageDesc *desc,
			 odouble doff[2]);

#endif /* OBITPOSITION_H  */
