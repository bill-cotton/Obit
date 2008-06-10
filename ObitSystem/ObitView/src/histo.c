/* $Id: histo.c,v 1.1 2005/07/23 00:42:27 bcotton Exp $ */
/* Histogram equalization and related functions  */ 
/*  Histogram equalization is an attempt to have equal numbers of pixels
    in each of the allowed color index states */
/*-----------------------------------------------------------------------
*  Copyright (C) 1998,1999-2008
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
  
#include <math.h>
#include "histo.h" 
#include "messagebox.h" 
#define WORKSIZE 4096 /* size of work histogram */

//*---------------Private function prototypes----------------*/
static olong find_hist (ObitFArray *pixels, ofloat lmax, ofloat lmin,
		ofloat **Hist);
    
/*---------------Public functions ----------------*/
/**
 * Get plane max and min pixel values
 * \param pixels array  for which to find extrema
 * \param max    maximum pixel value (blanked if failure)
 * \param min    minimum pixel value (blanked if failure)
 * \return returns 0 if OK else failed.
 */
olong get_extrema (ObitFArray *pixels, ofloat *max, ofloat *min) 
{
  olong pos[7];

  *max = ObitFArrayMax (pixels, &pos[0]);
  *min = ObitFArrayMin (pixels, &pos[0]);

  return 0;
} /* end get_extrema */

/**
 * Get plausible range of pixel values
 * \param pixels    array of pixel values
 * \param mapFunc   mapping type, 0=linear, 1=nonlinear, 
 *                  2=histogram equalization
 * \param max (input)  maximum value in image, if blanked determine
 *            (output) maximum pixel value to display
 * \param min (input)  minimum value in image, if blanked determine
 *            (output) minimum pixel value to display
 * \return returns 0 if OK else failed.
 */
olong get_range (ObitFArray *pixels, olong mapFunc,
		ofloat *max, ofloat *min) 
{
  ofloat lmax, lmin, lmode, lrms;
  ofloat blanked = ObitMagicF();

  if (!pixels) return -1; /* sanity check */
  if (!max) return -1;
  if (!min) return -1;

  /* get image range */
  lmax = *max;
  lmin = *min;
  if ((lmax==blanked) || (lmin==blanked) || (lmin>=lmax)) {
    get_extrema(pixels, &lmax, &lmin);
  }

  /* get rms, mode */
  lrms  = ObitFArrayRMS(pixels);
  lmode = ObitFArrayMode(pixels);
  
  /* set range by mapping type */
  if (mapFunc==0) { /* linear*/
    *min = lmode - lrms;
    *max = lmode + 0.1 * MAXCOLOR * lrms;
 } else if (mapFunc==1) { /* non linear */
    *min = lmode;
    *max = lmode + 0.3 * MAXCOLOR * lrms;
  } else if (mapFunc==2) { /* histogram equalization*/
    *min = lmode - lrms;
    *max = lmode + 0.5 * MAXCOLOR * lrms;
  }
  /* clip to range in image */
  if (*min<lmin) *min = lmin;
  if (*max>lmax) *max = lmax;

  return 0;
} /* end get_range */

/**
 * Compute histogram equalization mapping function, this must be called
 * before map_pixel.
 * \param pixels    array of pixel values
 * \param max (input)  maximum value in image, if blanked determine
 *            (output) maximum pixel value to display
 * \param min (input)  minimum value in image, if blanked determine
 *            (output) minimum pixel value to display
 * \param newRange TRUE iff a new display range to be determined
 * \param eq_map   Histogram equalization map for map_pixel, allocated
 *                 should be g_freeed when done.
 * \return returns 0 if OK else failed.
 */
int equalize (ObitFArray *pixels, float *max, float *min, ofloat **eq_map) 
{
  int   i, k, iHb, iHe;
  ofloat *hist=NULL;
  ofloat count, lmax, lmin, range, irange, *map;
  ofloat sum, perColor, fact, blanked = ObitMagicF();
  
  if (!pixels) return -1; /* sanity check */
  if (!max) return -1;
  if (!min) return -1;
  
  /* get clipping range */
  lmax = *max;
  lmin = *min;
  if ((lmax==blanked) || (lmin==blanked) || (lmin>=lmax)) {
    get_extrema(pixels, &lmax, &lmin);
  }

  count = (float) find_hist (pixels, lmax, lmin, &hist);
  /* if this failed it's probably still OK */
  if (count<=0.0) return -1;

  /* check for minimum number of pixels (10) */
  if (count<10.0) {
    /* Clean up and bail out on error */
    MessageShow ("Too few pixels in Pixel Range - Please reset");
    if (hist!=NULL) g_free(hist);
    return -3;
  } /* end of minimum count check */

  /* init histogram factors */
  range = lmax-lmin;
  if (range!=0.0) irange = (float)(WORKSIZE) / range;
  else irange =  (float)(WORKSIZE);
  fact = 1.0 / irange;

  iHb = 0;
  iHe = WORKSIZE-1;
  /* convert to mapping function */

  /* Allocate mapping function memory */
  map = g_malloc0((MAXCOLOR+1)*sizeof(ofloat));
  if (*eq_map) g_free(*eq_map);
  *eq_map = map;

  /* loop over function */
  perColor = count / ((ofloat)MAXCOLOR-1.0); /* pixels per color index */
  i = iHb;
  sum = 0.0;
  for (k=1; k<MAXCOLOR-1; k++) {
    while ((sum<perColor) && (i<=iHe)) {sum += hist[i++];}
    sum -= perColor; /* if more than allowed in one bin */
    *(map++) = lmin + (ofloat) i * fact;
  } /* end of loop over mapping function */

  /* first and last cells are the extrema */
  (*eq_map)[0]          = lmin;
  (*eq_map)[MAXCOLOR-1] = lmax;

  /* cleanup */
  if (hist) g_free(hist);
  return 0;
} /* end equalize */

/**
 * Return color index for specified pixel value.
 * 
 * \param map   histogram equalizatioin map computed by equalize
 * \param value Image pixel value
 * \return  color index , 0 used only for blanked values
 */
olong map_pixel(ofloat *map, float value) 
{
  int next, skip;
  
  if (!map) return 0; /* sanity check */
  
  /* clip to range */
  if (value<=map[0]) return 1;
  if (value>=map[MAXCOLOR-1]) return MAXCOLOR-1;
  
  /* lookup in table  - use binary search*/
  next = (MAXCOLOR / 2) - 1;
  skip = MAXCOLOR/4;
  while (skip>0) {
    if (value>map[next]) next += skip;
    else next -= skip;
    skip = skip / 2;
  }
  if (value<map[next]) next--;
  if (next<1) next = 1; /* 0 reserved for blank */
  return next;
  /* for (i=2; i<MAXCOLOR; i++) if (value < *map++) return i-1;     */
  
  return MAXCOLOR-1; /* just in case it gets here */
} /* end  map_pixel */

/*---------------Private functions ----------------*/

/**
 * Compute histogram clipping to [min,max]
 * \param pixels array  for which to find histogram
 * \param max    maximum allowed pixel value
 * \param min    minimum allowed pixel value
 * \param Hist   Histogram array, will be allocated, previous value
 *               will be deallocated.
 * \return returns 0 if OK else failed.-1=input error, 
 *         -2=allocation error 
 */
static olong find_hist (ObitFArray *pixels, float lmax, float lmin,
		 ofloat **Hist) 
{
  ofloat val, *hist, *valp, blanked = ObitMagicF();
  olong   i, j, k, nx, ny;
  olong   index;
  ofloat range, irange;
  olong  count;
  
  if (!pixels) return -1; /* sanity check */
  if (!Hist) return -1;
  
  /* allocate a work array for 4096 histogram levels */
  if (*Hist!=NULL) g_free(*Hist);
  *Hist = g_malloc0(WORKSIZE*sizeof(ofloat));
  if (*Hist==NULL) return -2;    /* allocation failure */
  hist = (ofloat*)*Hist;

  /* get dimensions */
  nx = pixels->naxis[0];
  ny = pixels->naxis[1];
  
  /* init histogram stuff */
  range = lmax-lmin;
  if (range!=0.0) irange = (ofloat)(WORKSIZE) / range;
  else irange =  (ofloat)(WORKSIZE);
  count = 0;
  for (k=0; k<WORKSIZE; k++) hist[k] = 0.0;

  /* loop over image */
  valp = pixels->array;  /* pointer in pixel array */
  for (j = 0; j<ny; j++) {
    for (i = 0; i<nx; i++) {
      val = *valp++;
      if (val!=blanked) {  /* only valid values */
	     index = (olong) ((irange * (val - lmin)) + 0.5);
	     if ((index>=0) && (index<WORKSIZE)) {
  	       count++;                /* total number of pixels */
	       hist[index] += 1.0;     /* accumulate histogram */
        }
      } /* end of valid pixel */
    } /* end loop over x */
  } /* end loop over y */

  return count;
} /* end find_hist */

