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
#include "Obit.h"
#ifndef OBITPOSLABELUTIL
#define OBITPOSLABELUTIL
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPosLabelUtil.h
 *  ObitPosLabelUtil Position labeling utility routine definition.
 */

/*---------------Public functions---------------------------*/

/** Public: appropriate character string describing a location */
void ObitPosLabelUtilAxisLabel(odouble pos, gchar* axis, gchar* label);

/** Public:  Convert RA in degrees to hh mm ss.sss */
void ObitPosLabelUtilRA2HMS(odouble ra, gchar* rach, gchar* rast);

/** Public: Convert dec in degrees to Dec dd mm ss.sss */
void ObitPosLabelUtilDec2DMS(odouble dec, gchar* decch, gchar* decst);

/** Public: Convert RA in degrees to hours min and seconds */
void ObitPosLabelUtilRAHMS (odouble ra, olong *h, olong *m, ofloat *s);

/** Public: Convert dec in degrees to degrees, min, sec */
void ObitPosLabelUtilDecDMS (odouble dec, int *d, olong *m, ofloat *s);

/** Public: Convert RA in hours min and seconds to degrees*/
olong ObitPosLabelUtilHMSRA (gint h, olong m, ofloat s, odouble *ra);

/** Public: convert dec in degrees, min, sec  to degrees */
olong ObitPosLabelUtilDMSDec (gint d, olong m, ofloat s, odouble *dec);
#endif /* OBITPOSLABELUTIL */ 
