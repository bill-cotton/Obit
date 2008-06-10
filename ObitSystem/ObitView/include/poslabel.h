/* $Id: poslabel.h,v 1.1 2005/07/23 00:20:08 bcotton Exp $ */
/*-----------------------------------------------------------------------
*  Copyright (C) 2005, 2008
*  Associated Universities, Inc. Washington DC, USA.
*  This program is free software; you can redistribute it and/or
*  modify it under the terms of the GNU General Public License as
*  published by the Free Software Foundation; either version 2 of
*  the License, or (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*-----------------------------------------------------------------------*/
#include "glib.h"
#include "ObitTypes.h"
#ifndef POSLABEL
#define POSLABEL
void AxisLabel(odouble pos, gchar* axis, gchar* label);
void ra2hms(odouble ra, gchar* rach, gchar* rast);
void dec2dms(odouble dec, gchar* decch, gchar* decst);
void rahms(odouble ra, olong *h, olong *m, ofloat *s);
void decdms(odouble dec, int *d, olong *m, ofloat *s);
int hmsra(gint h, olong m, ofloat s, odouble *ra);
int dmsdec(gint d, olong m, ofloat s, odouble *dec);
#endif /* POSLABEL */ 
