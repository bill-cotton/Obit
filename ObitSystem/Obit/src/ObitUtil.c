/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2008                                               */
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

#include "ObitUtil.h"
#include <gsl/gsl_multifit.h>

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUtil.c
 * ObitUtil utility function definitions.
 *
 */


/*----------------------Private functions---------------------------*/
/** Private: qsort ofloat comparison */
static int compare_gfloat  (const void* arg1,  const void* arg2);


/*----------------------Public functions---------------------------*/
/**
 * Get mean values of an array.
 * \param array array of values may be magic value blanked
 * \param n     dimension of array
 * \return      Mean value, possibly blanked
 */
ofloat meanValue (ofloat *array, olong incs, olong n)
{
  olong i, count;
  ofloat sum, fblank =ObitMagicF() ;

  if (n<=0) return fblank;
  sum = 0.0;
  count = 0;
    for (i=0; i<n; i++) {
      if (array[i*incs]!=fblank) { 
	sum += array[i*incs];
	count++;
      }
    }
    if (count>0) sum /= count;
    else sum = fblank;

  return sum;
} /* end meanValue */

/**
 * Get median value of an array.
 * Does sort then returns value array[ngood/2]
 * \param array   array of values, on return will be in ascending order
 * \param n       dimension of array
 * \param inc     stride in array
 * \return        Median value, possibly blanked
 */
ofloat medianValue (ofloat *array, olong incs, olong n)
{
  olong i, ngood;
  ofloat out, fblank = ObitMagicF();
  gboolean blanked;

  if (n<=0) return fblank;

  /* Count good (non blank) points */
  ngood = 0;
  for (i=0; i<n; i++) if (array[i*incs]!=fblank) ngood++;
  if (ngood<1) return fblank;  /* Any good data? */
  blanked = ngood<n;   /* Anything blanked? */

  /* set blanked values to 2.0e20 for sort */
  for (i=0; i<n; i++) if (array[i*incs]==fblank) array[i*incs] = 2.0e20;

  /* Sort to ascending order */
  qsort ((void*)array, n, MAX(1,incs)*sizeof(ofloat), compare_gfloat);

  /* Set median */
  out = array[(ngood/2)*incs];

  /* reset blanked values */
  if (blanked) {
    for (i=0; i<n; i++) if (array[i*incs]>1.0e20) array[i*incs] = fblank;
  }

  return out;
} /* end medianValue */

/**
 * Return average of navg values around median (with magic value blanking)
 * Does sort then returns average of center of array
 * If there are fewer than navg valid data then the median is returned.
 * \param array array of values, on return will be in ascending order
 * \param incs  increment in data array, >1 => weights
 * \param navg  width of averaging about median
 * \param doWt  If True value after data is a weight
 * \param n     dimension of array
 * \return      Median/average, possibly blanked
 */
ofloat medianAvg (ofloat *array, olong incs, olong navg, gboolean doWt, olong n)
{
  olong i, cnt, wid, ngood, ind, hi, lo;
  ofloat temp, wt, out=0.0, fblank = ObitMagicF();
  ofloat center;
  gboolean blanked;

  if (n<=0) return fblank;

  /* Count good (non blank) points */
  ngood = 0;
  for (i=0; i<n; i++) if (array[i*incs]!=fblank) ngood++;
  if (ngood<1) return fblank;  /* Any good data? */
  blanked = ngood<n;   /* Anything blanked? */

  /* set blanked values to 2.0e20 for sort */
  for (i=0; i<n; i++) if (array[i*incs]==fblank) array[i*incs] = 2.0e20;

  /* Sort to ascending order */
  qsort ((void*)array, n, MAX(1,incs)*sizeof(ofloat), compare_gfloat);

  /* How many good values? */
  cnt = 0;
  for (i=0; i<n; i++) if (array[i*incs]<1.0e20) cnt++;
  wid = MIN (MAX (1, navg/2), ngood/2);  /* width of averaging */
  ngood = cnt;

  ind = (ngood-1)/2;  /* center good value */
  out = array[ind*incs];

  /* average central values */
  if (ngood>navg) {
    center = (ngood-1)/2.0;
    lo = MAX (1, (olong)(center - wid + 0.6));
    hi = MIN (ngood, (olong)(center + wid));
    /* Weighted? */
    if (doWt) {
      temp = 0.0; wt = 0.0;
      for (i=lo; i<=hi; i++) {
	if (array[i*incs]<1.0e20) {
	  wt += array[i*incs+1];
	  temp += array[i*incs]*array[i*incs+1];
	}
      }
      if (wt>0.0) out = temp/wt;
    } else { /* unweighted */
      temp = 0.0; cnt = 0;
      for (i=lo; i<=hi; i++) {
	if (array[i]<1.0e20) {
	  cnt++;
	  temp += array[i];
	}
	if (cnt>0) out = temp/cnt;
      }
    } /* end if weighting */
  }

  /* reset blanked values */
  if (blanked) {
    for (i=0; i<n; i++) if (array[i*incs]>1.0e20) array[i*incs] = fblank;
  }

  return out;
} /* end medianAvg */

/**
 * Return the running median of an array
 * \param n       Number of points
 * \param wind    Width of median window in cells
 * \param array   Array of values, fblank blanking supported
 * \param alpha   0 -> 1 = pure boxcar -> pure MWF (Alpha of the 
 *                data samples in a window are discarded and 
 *                the rest averaged). 
 * \param rms     RMS of array, median average sigma from each wind of data
 * \param out     array of size of array to be with median values
 * \param work    work array of size of array
 */
void RunningMedian (olong n, olong wind, ofloat *array, ofloat alpha, 
		    ofloat *RMS, ofloat *out, ofloat *work)
{
  ofloat *lwork=NULL;
  ofloat level, sigma, sigmaSum, sigmaCnt;
  ofloat fblank = ObitMagicF();
  olong i, j, k, op, ind, half, RMScnt=0;

  /* Create array */
  lwork = g_malloc0(wind*sizeof(ofloat));

  half = wind/2;
  ind  = 0;
  op   = 0;
  k    = 0;
  sigmaSum = 0.0;
  sigmaCnt = 1;

  /* First half wind filled with median of first wind points */
  for (j=ind; j<ind+wind; j++) lwork[k++] = array[j];
  ind++;
  level = MedianLevel (wind, lwork, alpha);
  sigma = MedianSigma (wind, lwork, level);
  if (sigma!=fblank) {
    sigmaSum = sigma;
    sigmaSum += sigma;
    work[RMScnt++] = sigma;
  }
  for (k=0; k<half; k++) out[op++] = level;

  /* Loop over middle of array */
  for (i=half; i<n-half; i++) {
    k = 0;
    for (j=ind; j<ind+wind; j++) lwork[k++] = array[j];
    ind++;
    level = MedianLevel (wind, lwork, alpha);
    sigma = MedianSigma (wind, lwork, level);
    if (sigma!=fblank) {
      sigmaSum += sigma;
      sigmaCnt++;
      work[RMScnt++] = sigma;
    }
    out[op++] = level;
  } /* end loop over array */

  /* Fill in bit at end */
  while (op<n) {
    out[op++] = level;
  }

  /* median Average sigmas */
  *RMS = MedianLevel (RMScnt, work, alpha);

  /* Cleanup */
  if (lwork) g_free(lwork);

} /* end RunningMedian */

/**
 * Determine alpha median value of a ofloat array
 * \param n       Number of points
 * \param value   Array of values
 * \param alpha   0 -> 1 = pure boxcar -> pure MWF (ALPHA of the 
 *                data samples are discarded and the rest averaged). 
 * \return alpha median value
 */
ofloat MedianLevel (olong n, ofloat *value, ofloat alpha)
{
  ofloat out=0.0;
  ofloat fblank = ObitMagicF();
  ofloat beta, sum;
  olong i, i1, i2, count;

  if (n<=0) return out;

  /* Sort to ascending order */
  qsort ((void*)value, n, sizeof(ofloat), compare_gfloat);

  out = value[n/2];

  beta = MAX (0.05, MIN (0.95, alpha)) / 2.0; /*  Average around median factor */

  /* Average around the center */
  i1 = MAX (0, (n/2)-(olong)(beta*n+0.5));
  i2 = MIN (n, (n/2)+(olong)(beta*n+0.5));

  if (i2>i1) {
    sum = 0.0;
    count = 0;
    for (i=i1; i<i2; i++) {
      if (value[i]!=fblank) {
	sum += value[i];
	count++;
      }
    }
    if (count>0) out = sum / count;
  }
   
  return out;
} /* end MedianLevel */

/**
 * Determine robust RMS value of a ofloat array about mean
 * Use center 90% of points, excluding at least one point from each end
 * \param n       Number of points, needs at least 4
 * \param value   Array of values assumed sorted
 * \param mean    Mean value of value
 * \return RMS value, fblank if cannot determine
 */
ofloat MedianSigma (olong n, ofloat *value, ofloat mean)
{
  ofloat fblank = ObitMagicF();
  ofloat out;
  ofloat sum;
  olong i, i1, i2, count;

  out = fblank;
  if (n<=4) return out;
  if (mean==fblank) return out;

  /* Get RMS around the center 90% */
  i1 = MAX (1,   (n/2)-(olong)(0.45*n+0.5));
  i2 = MIN (n-1, (n/2)+(olong)(0.45*n+0.5));

  if (i2>i1) {
    sum = 0.0;
    count = 0;
    for (i=i1; i<i2; i++) {
      if (value[i]!=fblank) {
	sum += (value[i]-mean)*(value[i]-mean);
	count++;
      }
    }
    if (count>1) out = sqrt(sum / (count-1));
  }
   
  return out;
} /* end MedianSigma */

/**
 * Fit polynomial y = f(poly, x) with magic value blanking
 * Use gsl package.
 * \param poly   [out] polynomial coef in order of increasing power of x
 * \param order  order of the polynomial
 * \param x      values at which y is sampled
 * \param y      values to be fitted
 * \param wt     weights for values
 * \param n      number of (x,y) value pairs
 */
void  FitPoly (ofloat *poly, olong order, ofloat *x, ofloat *y, ofloat *wt, 
	       olong n)
{
  olong i, j, k, good, p=order+1;
  ofloat fgood=0.0;
  double xi, chisq;
  gsl_matrix *X, *cov;
  gsl_vector *yy, *w, *c;
  gsl_multifit_linear_workspace *work;
  ofloat fblank = ObitMagicF();

  /* Only use good data */
  good = 0;
  for (i=0; i<n; i++) if (y[i]!=fblank) {good++; fgood = y[i];}

  /* If only one good datum - use it */
  if (good==1) {
    poly[0] = fgood;
    poly[1] = 0.0;
  }

  if (good<(order+1)) {  /* Not enough good data */
    poly[0] = fblank;
    return;
  }

  /* If only one and order=0 use it */
  if ((good==1) && (order==0)) {
    poly[0] = y[0];
    return;
  }

  /* allocate arrays */
  X    = gsl_matrix_alloc(good, p);
  yy   = gsl_vector_alloc(good);
  w    = gsl_vector_alloc(good);
  c    = gsl_vector_alloc(p);
  cov  = gsl_matrix_alloc(p, p);
  work = gsl_multifit_linear_alloc (good, p);

  /* set data */
  k = 0;
  for (i=0; i<n; i++) {
    if (y[i]!=fblank) {
      gsl_vector_set(yy, k, (double)y[i]);
      gsl_vector_set(w, k, (double)wt[i]);
      xi = 1.0;
      for (j=0; j<p; j++) {
	gsl_matrix_set(X, k, j, xi);
	xi *= x[i];
      }
      k++;  /* good ones */
    }
  }

  /* Fit */
  gsl_multifit_wlinear (X, w, yy, c, cov, &chisq, work);

  /* get results */
  for (j=0; j<p; j++) poly[j] = (ofloat)gsl_vector_get(c, j);

  /* Deallocate arrays */
  gsl_matrix_free(X);
  gsl_vector_free(yy);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free (work);
} /* end FitPoly */

/**
 * Evaluate polynomial y = f(poly, x)
 * \param order  order of the polynomial
 * \param poly   polynomial coef in order of increasing power of x
 * \param x      value at which polynomial is to be evaluated
 * \return evaluated polynomial
 */
ofloat  EvalPoly (olong order, ofloat *poly, ofloat x)
{
  olong i;
  ofloat sum = 0.0, aarg = 1.0;

  for (i=0; i<=order; i++) {
    sum  += poly[i]*aarg;
    aarg *= x;
  }
  return sum;
} /* end EvalPoly  */



/*----------------------Private functions---------------------------*/
/**
 * ofloat comparison of two arguments
 * \param arg1 first value to compare
 * \param arg2 second value to compare
 * \return negative if arg1 is less than arg2, zero if equal
 *  and positive if arg1 is greater than arg2.
 */
static int compare_gfloat  (const void* arg1,  const void* arg2)
{
  int out = 0;
  ofloat larg1, larg2;

  larg1 = *(ofloat*)arg1;
  larg2 = *(ofloat*)arg2;
  if (larg1<larg2) out = -1;
  else if (larg1>larg2) out = 1;
  return out;
} /* end compare_gfloat */
