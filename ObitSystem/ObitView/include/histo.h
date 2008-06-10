/* $Id: histo.h,v 1.1 2005/07/23 00:42:27 bcotton Exp $ */
/*  header file for histogram equalization functions for Matrix Class */ 
/*  Histogram equalization is an attempt to have equal numbers of pixels
    in each of the allowed color index states */
/*-----------------------------------------------------------------------
*  Copyright (C) 1998-2008
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
  
#include <glib.h> 
#include "ObitFArray.h" 
#ifndef HISTO_H
#define HISTO_H 
#define MAXCOLOR 128     /* number of colors in display 
			  also defined in obitview.h */
/* get plane max and min pixel values
   pixels = array for which to find extrema
   max (returned) = maximum pixel value.
   min (returned) = minimum pixel value
   Returns 0 if worked else failed */
olong get_extrema (ObitFArray *pixels, ofloat *max, ofloat *min);

/* get plausible range of pixel values
   pixels = array for which to find extrema
   max (returned) = maximum pixel value.
   min (returned) = minimum pixel value
   Returns 0 if worked else failed */
olong get_range (ObitFArray *pixels, olong mapFunc,
		ofloat *max, ofloat *min);

/* compute histogram equalization mapping function, this must be called
   before map_pixel;  after the image is opened (LoadFImage) and before
   the first read (Readpatch).
   pixels = array for which to make equalize histogram
   max = maximum pixel value (clip above)
   min = minimum pixel value (clip below)
   eq_map = histogram equalization map for map_pixel
   Returns 0 if worked else failed */
olong equalize (ObitFArray *pixels, ofloat *max, ofloat *min, ofloat **eq_map);

/* Return color index for specified pixel value. Returns 0 (blanked) if 
   mapping function invalid.
   pixels = array for which to make equalize histogram,
   value = pixel value,
   returns color index (0 used only for blanked values */
olong map_pixel(ofloat *eq_map, ofloat value);
  
#endif /* HISTO_H */ 
