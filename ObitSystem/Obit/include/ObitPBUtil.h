/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2014                                          */
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
#ifndef OBITPBUTIL_H 
#define OBITPBUTIL_H 

#include "Obit.h"
#include "ObitImage.h"
#include "ObitTableCC.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPBUtil.h
 * Antenna primary beam shape utility module.
 *
 */


/*---------------Public functions---------------------------*/
/** Use polynomial beam shape - useful for VLA frequencies < 1.0 GHz */
ofloat ObitPBUtilPoly (odouble Angle, odouble Freq, ofloat pbmin);

/** Use Jinc beam shape - useful for frequencies > 1.0 GHz */
ofloat ObitPBUtilJinc (odouble Angle, odouble Freq, ofloat antSize, ofloat pbmin);

/** Use KAT-7 beam shape  */
ofloat ObitPBUtilKAT7 (odouble Angle, odouble Freq, ofloat pbmin);

/** Function which returns relative primary beam correction */
ofloat ObitPBUtilRelPB (odouble Angle, olong nfreq, odouble *Freq, ofloat antSize, 
			ofloat pbmin, odouble refFreq);

/** Function which returns pointing error amplitude correction */
ofloat ObitPBUtilPntErr (odouble Angle, odouble AngleO, ofloat antSize, 
			 ofloat pbmin, odouble Freq);

/** Correct ObitTableCC for relative Primary Beam */
ObitTableCC *ObitPBUtilCCCor(ObitImage *image, olong inCCver, olong *outCCver, 
			     olong nfreq, odouble *Freq, ofloat antSize, ofloat pbmin,
			     odouble refFreq, olong *startCC, olong *endCC, 
			     ObitErr *err);

/** Correct Image for relative Primary Beam */
ObitFArray* ObitPBUtilImageCor(ObitImage *image, olong *inPlane, 
			       olong nfreq, odouble *Freq, 
			       ofloat antSize, ofloat pbmin, odouble refFreq, 
			       ObitErr *err);

#endif /* OBITPBUTIL_H */ 
