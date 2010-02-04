/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2010                                          */
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

#if HAVE_GSL==1  /* GSL stuff */
#include "gsl/gsl_blas.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_multifit_nlin.h"
#endif /* GSL stuff */
#include "ObitFArrayUtil.h"
#include "ObitCArray.h"
#include "ObitFFT.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFArrayUtil.c
 * ObitFArray utility function definitions.
 */

/*---------------Private function prototypes----------------*/
/** Private: Fit 2-d circular Gaussian */
static ofloat FitCGauss (olong Count, ofloat *pixX, ofloat *pixY, ofloat *val, 
			 ofloat *peak, ofloat *center, ofloat *sigma, 
			 ObitErr *err);

/*----------------------Public functions---------------------------*/
/**
 * Fits 2-D circular Gaussian in an FArray.
 * Starting solution is peak value in array in.
 * \param in      2-D FArray to fit
 * \param FWHM    [in]  half width of box around peak to use in fitting
 *                      0 => all.
 * \param FWHM    [out] Full width malf max of Gaussian
 * \param center  [out] 0-rel center pixel position
 * \param peak    [out] Peak value
 * \param err      Error stack, returns if not empty.
 * \return RMS residual, -1 on error
 */
ofloat ObitFArrayUtilFitCGauss (ObitFArray *in, ofloat *FWHM, ofloat *center, 
				ofloat *peak, ObitErr *err)
{
  ofloat RMS = -1.0;
  ofloat *pixX=NULL, *pixY=NULL, *val=NULL, sigma;
  ofloat *data, fblank = ObitMagicF();
  olong size, Count, ix, iy, nx, ny, pos[2];
  olong blc[2]= {0,0}, trc[2]={0,0}, ihalf;
  gboolean subBox;
  ObitFArray *box = NULL;
  gchar *routine = "ObitFArrayUtilFitCGauss";

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return RMS;
  g_assert (ObitFArrayIsA(in));

  /* Need subarray on in? */
  ihalf = (olong)((*FWHM) + 0.5);
  subBox = ihalf > 0;

  /* Initial guess from peak */
  *peak = ObitFArrayMax (in, pos);
  center[0] = (ofloat)pos[0];
  center[1] = (ofloat)pos[1];
  *FWHM = 4.0;

  if (subBox) {
    blc[0] = MAX (0, pos[0]-ihalf);
    blc[1] = MAX (0, pos[1]-ihalf);
    trc[0] = MIN (in->naxis[0]-1, pos[0]+ihalf);
    trc[1] = MIN (in->naxis[1]-1, pos[1]+ihalf);
    box = ObitFArraySubArr (in, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, RMS);
    *peak = ObitFArrayMax (box, pos);
    center[0] = (ofloat)pos[0];
    center[1] = (ofloat)pos[1];
    nx    = box->naxis[0];
    ny    = box->naxis[1];
    data  = box->array;
  } else { /* use all */
    nx    = in->naxis[0];
    ny    = in->naxis[1];
    data  = in->array;
  }

  /* Create arrays */
  size = nx * ny;
  pixX = g_malloc0(size*sizeof(ofloat));
  pixY = g_malloc0(size*sizeof(ofloat));
  val  = g_malloc0(size*sizeof(ofloat));
  
  /* Data from array */
  Count = 0;
  for (iy=0; iy<ny; iy++) {
    for (ix=0; ix<nx; ix++) {
      if (data[ix+iy*nx] != fblank ) {
	pixX[Count] = (ofloat)ix;
	pixY[Count] = (ofloat)iy;
	val[Count]  = data[ix+iy*nx];
	Count++;
      }
    }
  }

  /* Fit */
  sigma = (*FWHM) / 2.355;  /* Use Gaussian sigma in fitting */
  RMS = FitCGauss (Count, pixX, pixY, val, peak, center, &sigma, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, RMS);
  *FWHM = 2.355 * sigma;

  /* Correct center for windowing */
  center[0] += blc[0];
  center[1] += blc[1];

  /* Cleanup */
  if (pixX) g_free(pixX);
  if (pixY) g_free(pixY);
  if (val)  g_free(val);
  if (box)  box = ObitFArrayUnref(box);
 
  return RMS;
}  /* end ObitFArrayUtilFitCGauss */

/**
 * Convolves two 2-D Arrays using FFTs
 * Arrays must have the same geometry and NOT contain magic value blanking
 * \param in1  First input array
 * \param in2  Second input array
 * \param err  Error stack, returns if not empty.
 * \return Convolution, shoud be ObitFArrayUnrefed when done
 */
ObitFArray* ObitFArrayUtilConvolve (ObitFArray *in1, ObitFArray *in2, 
				    ObitErr *err)
{
  ObitFArray *out=NULL;
  ObitCArray *uv1 = NULL, *uv2 = NULL;
  ObitFFT *forFFT = NULL, *revFFT = NULL;
  olong  ddim[2];
  olong ndim, naxis[2];
  ofloat scale;
  gchar *routine = "ObitFArrayUtilConvolve";

   /* error checks */
  if (err->error) return out;
  g_assert (ObitFArrayIsA(in1));
  g_assert (ObitFArrayIsA(in2));
  Obit_retval_if_fail(ObitFArrayIsCompatable(in1, in2),  err, out,
		      "%s: FArrays %s and %s incompatible", 
		      routine, in1->name, in2->name);
  Obit_retval_if_fail((in1->ndim==2),  err, out,
		      "%s: FArrays %s and %s NOT 2-D", 
		      routine, in1->name, in2->name);
  Obit_retval_if_fail((in1->naxis[0]==ObitFFTSuggestSize(in1->naxis[0])),  
		      err, out,
		      "%s: FArray %s dim 1 NOT proper size for FFT", 
		      routine, in1->name);
  Obit_retval_if_fail((in1->naxis[1]==ObitFFTSuggestSize(in1->naxis[1])),  
		      err, out,
		      "%s: FArray %s dim 2NOT proper size for FFT", 
		      routine, in1->name);

  /* Make UV plane arrays */
  ndim = 2;
  naxis[0] = 1+in1->naxis[0]/2; naxis[1] = in1->naxis[1]; 
  uv1 = ObitCArrayCreate ("Convolve work 1", ndim, naxis);
  uv2 = ObitCArrayCreate ("Convolve work 2", ndim, naxis);

  /* Make FFTs */
  ddim[0] = in1->naxis[0]; ddim[1] = in1->naxis[1];
  forFFT = newObitFFT("FTImage", OBIT_FFT_Forward, 
		      OBIT_FFT_HalfComplex, 2, ddim);
  revFFT = newObitFFT("FTuv", OBIT_FFT_Reverse, 
		      OBIT_FFT_HalfComplex, 2, ddim);

  /* FFT to uv plane */
  ObitFArray2DCenter (in1); /* Swaparoonie to FFT order */
  ObitFFTR2C (forFFT, in1, uv1);
  ObitFArray2DCenter (in2);
  ObitFFTR2C (forFFT, in2, uv2);

  /* return input to original order */
  ObitFArray2DCenter (in1);
  ObitFArray2DCenter (in2);

  /* Scale */
  scale = 1.0 / ((ofloat)in1->naxis[0] * (ofloat)in1->naxis[1]);
  ObitCArraySMul (uv1, scale);
  /*ObitCArraySMul (uv2, scale); only 1?*/

  /* Multiply */
  ObitCArrayMul(uv1, uv2, uv1);

  /* Some Cleanup */
  forFFT = ObitFFTUnref(forFFT);
  uv2    = ObitCArrayUnref(uv2);

  /* Create output */
  out =  newObitFArray ("Convolution");;
  ObitFArrayClone (in1, out, err);
  if (err->error) {
    uv1    = ObitCArrayUnref(uv1); /* cleanup */
    revFFT = ObitFFTUnref(revFFT);
    Obit_traceback_val (err, routine, in1->name, out);
  }

  /* FFT back to image plane */
  ObitFFTC2R (revFFT, uv1, out);
  
  /* Put the center at the center */
  ObitFArray2DCenter (out);

  uv1    = ObitCArrayUnref(uv1); /* cleanup */
  revFFT = ObitFFTUnref(revFFT);

  return out;
} /* end ObitFArrayUtilConvolve */

/**
 * Correlate two 2-D Arrays using FFTs
 * Arrays must have the same geometry and NOT contain magic value blanking
 * \param in1  First input array
 * \param in2  Second input array
 * \param err  Error stack, returns if not empty.
 * \return Convolution, shoud be ObitFArrayUnrefed when done
 */
ObitFArray* ObitFArrayUtilCorrel (ObitFArray *in1, ObitFArray *in2, 
				  ObitErr *err)
{
  ObitFArray *out=NULL;
  ObitCArray *uv1 = NULL, *uv2 = NULL;
  ObitFFT *forFFT = NULL, *revFFT = NULL;
  olong  ddim[2];
  olong ndim, naxis[2];
  ofloat scale;
  gchar *routine = "ObitFArrayUtilCorrel";

   /* error checks */
  if (err->error) return out;
  g_assert (ObitFArrayIsA(in1));
  g_assert (ObitFArrayIsA(in2));
  Obit_retval_if_fail(ObitFArrayIsCompatable(in1, in2),  err, out,
		      "%s: FArrays %s and %s incompatible", 
		      routine, in1->name, in2->name);
  Obit_retval_if_fail((in1->ndim==2),  err, out,
		      "%s: FArrays %s and %s NOT 2-D", 
		      routine, in1->name, in2->name);
  Obit_retval_if_fail((in1->naxis[0]==ObitFFTSuggestSize(in1->naxis[0])),  
		      err, out,
		      "%s: FArray %s dim 1 NOT proper size for FFT", 
		      routine, in1->name);
  Obit_retval_if_fail((in1->naxis[1]==ObitFFTSuggestSize(in1->naxis[1])),  
		      err, out,
		      "%s: FArray %s dim 2NOT proper size for FFT", 
		      routine, in1->name);

  /* Make UV plane arrays */
  ndim = 2;
  naxis[0] = 1+in1->naxis[0]/2; naxis[1] = in1->naxis[1]; 
  uv1 = ObitCArrayCreate ("Correl work 1", ndim, naxis);
  uv2 = ObitCArrayCreate ("Correl work 2", ndim, naxis);

  /* Make FFTs */
  ddim[0] = in1->naxis[0]; ddim[1] = in1->naxis[1];
  forFFT = newObitFFT("FTImage", OBIT_FFT_Forward, 
		      OBIT_FFT_HalfComplex, 2, ddim);
  revFFT = newObitFFT("FTuv", OBIT_FFT_Reverse, 
		      OBIT_FFT_HalfComplex, 2, ddim);

  /* FFT to uv plane */
  ObitFArray2DCenter (in1); /* Swaparoonie to FFT order */
  ObitFFTR2C (forFFT, in1, uv1);
  ObitFArray2DCenter (in2);
  ObitFFTR2C (forFFT, in2, uv2);

  /* return input to original order */
  ObitFArray2DCenter (in1);
  ObitFArray2DCenter (in2);

  /* Scale */
  scale = 1.0 / ((ofloat)in1->naxis[0] * (ofloat)in1->naxis[1]);
  ObitCArraySMul (uv1, scale);
  /*ObitCArraySMul (uv2, scale); only 1?*/

  /* conjugate uv2 */
  ObitCArrayConjg(uv2);

  /* Multiply */
  ObitCArrayMul(uv1, uv2, uv1);

  /* Some Cleanup */
  forFFT = ObitFFTUnref(forFFT);
  uv2    = ObitCArrayUnref(uv2);

  /* Create output */
  out =  newObitFArray ("Convolution");;
  ObitFArrayClone (in1, out, err);
  if (err->error) {
    uv1    = ObitCArrayUnref(uv1); /* cleanup */
    revFFT = ObitFFTUnref(revFFT);
    Obit_traceback_val (err, routine, in1->name, out);
  }

  /* FFT back to image plane */
  ObitFFTC2R (revFFT, uv1, out);
  
  /* Put the center at the center */
  ObitFArray2DCenter (out);

  uv1    = ObitCArrayUnref(uv1); /* cleanup */
  revFFT = ObitFFTUnref(revFFT);

  return out;
} /* end ObitFArrayUtilCorrel */

/**
 * Create a Gaussian UV tapering array corresponding to an image plane Gaussian
 *  Lifted from AIPS/CONVL.FOR
 * \param naxis   Dimension of image as [nx,ny]
 * \param cells   Cell spacing in x and y in units of maj,min (asec)
 * \param maprot  Map rotation (deg)
 * \param maj     Major axis of Gaussian in image plane (same units as cells)
 * \param min     Minor axis of Gaussian in image plane (same units as cells)
 * \param pa      Position angle of Gaussian in image plane, from N thru E, (deg)
 * \return ObitFArray, should be unReffed when done, as [u,v]
 */
ObitFArray* ObitFArrayUtilUVGaus (olong *naxis, ofloat *cells, ofloat maprot,
				  ofloat Gaumaj, ofloat Gaumin, ofloat GauPA)
{
  ObitFArray* outFA = NULL;
  ofloat ta, tb, am, an, xnx2, xny2, xnxny, gausaa, gausbb, gauscc;
  ofloat arg, varg, amp, u, v;
  olong iu, iv, ucen, vcen;

  /* Create output */
  outFA = ObitFArrayCreate("UVpix",2, naxis);

  ta = Gaumaj * G_PI / 1.1774;
  tb = Gaumin * G_PI / 1.1774;
  am = cos((GauPA+maprot-90.)*DG2RAD);
  an = sin((GauPA+maprot-90.)*DG2RAD);
  xnx2 = naxis[0] * cells[0];
  xny2 = naxis[1] * cells[1];
  xnxny = fabsf (xnx2 * xny2);
  xnx2 *= xnx2;
  xny2 *= xny2;
  gausaa = 0.5*(ta*ta*am*am + tb*tb*an*an) / (xny2);
  gausbb = ((tb*tb-ta*ta) * an*am) / (xnxny);
  gauscc = 0.5*(ta*ta*an*an + tb*tb*am*am) / (xnx2);

  ucen = naxis[0]/2;
  vcen = naxis[1]/2;
  /* Loop computing Gaussian */
  for (iv=0; iv<naxis[1]; iv++) {
    v = (iv - vcen);
    varg = v*v*gauscc;
    for (iu=0; iu<naxis[0]; iu++) {
      u = (iu - ucen);
      arg = -(u*u*gausaa + u*v*gausbb + varg);
      amp = expf(arg);
      outFA->array[iv*naxis[0]+iu] = amp;
    }
  }
  return outFA;
} /* end ObitFArrayUtilUVGaus */

/*---------------Private functions--------------------------*/

/* Structure for least squares fitting */
struct fitData {
  size_t n;   /* number of data points */
  float *x;   /* 1st independent variable */
  float *y;   /* 2nd independent variable */
  float *val; /* dependent variable */
  float RMS;  /* RMS residual */
};

#if HAVE_GSL==1  /* GSL stuff */
/**
 * CircularGaussian model calculating routine for gsl least squares fitter
 * \param coef   Coefficient array (amplitude, center[2], 1/variance)
 * \param params Data structure
 * \param f      [out] function residuals
 * \returns GSL completion code
 */
static int cgaussFunc (const gsl_vector *coef, void *params, gsl_vector *f)
{
  size_t n   = ((struct fitData *)params)->n;
  float *x   = ((struct fitData *)params)->x;
  float *y   = ((struct fitData *)params)->y;
  float *val = ((struct fitData *)params)->val;
  float sum = 0.0;
  long i;
  double amp, cenx, ceny, ivar, model, resid;

  /* Current parameters */
  amp   = gsl_vector_get(coef, 0);
  cenx  = gsl_vector_get(coef, 1);
  ceny  = gsl_vector_get(coef, 2);
  ivar  = gsl_vector_get(coef, 3);

  /* Loop through data calculating residuals to model */
  sum = 0.0;
  for (i=0; i<n; i++) {
    model = amp * 
      exp (-((x[i]-cenx)*(x[i]-cenx) + (y[i]-ceny)*(y[i]-ceny)) * ivar);
    resid = model - val[i];
    gsl_vector_set(f, i, resid);  /* to output vector */
    sum += resid*resid;
  }

  ((struct fitData *)params)->RMS = sqrt (sum/n);

  return GSL_SUCCESS;
} /* end  cgaussFunc */

/**
 * Circular Gaussian model calculating Jacobean for gsl least squares fitter
 * This is the partial derivative of the residuals matrix
 * \param coef   Coefficient array (amplitude, center[2], 1/variance)
 * \param params Data structure
 * \param J      [out] Jacobean values
 * \returns GSL completion code
 */
static int cgaussJacob (const gsl_vector *coef, void *params, gsl_matrix *J)
{
  size_t n   = ((struct fitData *)params)->n;
  float *x   = ((struct fitData *)params)->x;
  float *y   = ((struct fitData *)params)->y;
  long i;
  double amp, cenx, ceny, ivar, eterm;
  double part1, part2, part3, part4;

  /* Current parameters */
  amp   = gsl_vector_get(coef, 0);
  cenx  = gsl_vector_get(coef, 1);
  ceny  = gsl_vector_get(coef, 2);
  ivar  = gsl_vector_get(coef, 3);

  /* Loop through data calculating partial derivatives of residuals */
  for (i=0; i<n; i++) {
    eterm = exp (-((x[i]-cenx)*(x[i]-cenx) + (y[i]-ceny)*(y[i]-ceny)) * ivar);
    part1 = eterm;                            /* partial wrt amplitude */
    part2 = +amp*eterm*ivar*2.0*(x[i]-cenx);  /* partial wrt center x */
    part3 = +amp*eterm*ivar*2.0*(y[i]-ceny);  /* partial wrt center y */
    /* partial wrt 1/var  */
    part4 = -amp*eterm*((x[i]-cenx)*(x[i]-cenx) + (x[i]-cenx)*(x[i]-cenx));
    gsl_matrix_set(J, i, 0, part1);  /* to output matrix */
    gsl_matrix_set(J, i, 1, part2);
    gsl_matrix_set(J, i, 2, part3);
    gsl_matrix_set(J, i, 3, part4);
  }

  return GSL_SUCCESS;
} /* end  cgaussJacob */


/**
 * Compute both function and derivatives for  gsl least squares fitter
 * \param coef   Coefficient array (amplitude, center[2], sigma)
 * \param params Data structure
 * \param f      [out] function residuals
 * \param J      [out] Jacobean
 * \returns GSL completion code
 */
static int cgaussFuncJacob (const gsl_vector *coef, void *params, gsl_vector *f, 
			   gsl_matrix *J)
{
  cgaussFunc (coef, params, f);
  cgaussJacob(coef, params, J);
  return GSL_SUCCESS;
} /* end  cgaussFuncJacob */
#endif /* GSL stuff */

/**
 * Fit 2D circular Gaussian
 * \param Count  Number of samples to be fitted
 * \param pixX   X pixel coordinate
 * \param pixY   Y pixel coordinate
 * \param val    Measured values
 * \param peak   [in/out] peak value
 * \param center [in/out] center pixel of Gaussian [2]
 * \param sigma  [in/out] Gaussian sigma (width) in pixels
 * \param err      Error stack, returns if not empty.
 * \retyrn RMS residual (-1.0 on failure)
 */
static ofloat FitCGauss (olong Count, ofloat *pixX, ofloat *pixY, ofloat *val, 
			 ofloat *peak, ofloat *center, ofloat *sigma, 
			 ObitErr *err)
{
#if HAVE_GSL==1  /* GSL stuff */
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver* solver;
  gsl_multifit_function_fdf func;
  gsl_vector_view coef;
  struct fitData data;
  int status, iter;
  double coef_init[4];
  gchar *routine = "FitCGauss";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return -1.0;
  if (Count<=0) return -1.0;
  
  /* Fitter data structure */
  data.n   = Count;
  data.x   = pixX;;
  data.y   = pixY;
  data.val = val;

  /* initial guess */
  coef_init[0] = *peak;
  coef_init[1] = center[0];
  coef_init[2] = center[1];
  coef_init[3] = 1.0 / (2.0 * (*sigma) * (*sigma)); /* fit as 1/2*var */

  coef = gsl_vector_view_array (coef_init, 4);

  /* Create /fill function structure */
  func.f   = &cgaussFunc;      /* Compute function */
  func.df  = &cgaussJacob;     /* Compute Jacobian (derivative matrix) */
  func.fdf = &cgaussFuncJacob; /* Compute both function and derivatives */
  func.n = Count;              /* Number of data points */
  func.p = 4;                  /* number of parameters */
  func.params = &data;         /* Data structure */

  T = gsl_multifit_fdfsolver_lmsder;
  solver = gsl_multifit_fdfsolver_alloc(T, Count, 4);
  gsl_multifit_fdfsolver_set(solver, &func, &coef.vector);

  /* ready to rumble */
  iter = 0;
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);
    if ((status!=GSL_CONTINUE) && (status!=GSL_SUCCESS)) {/* problem? */
      Obit_log_error(err, OBIT_Error, "%s: Solver status %d %s", 
		     routine,status,gsl_strerror(status));
      break; 
    }
    /* convergence test */
    status = gsl_multifit_test_delta(solver->dx, solver->x, 1.0e-4, 1.0e-4);
   }
  while ((status == GSL_CONTINUE) && (iter< 500));

  /* return results */
  /* debug fprintf (stderr,"no. iter %d base %f peak %f\n",iter, *base,*peak);*/
  
  *peak     = (float)gsl_vector_get(solver->x, 0);
  center[0] = (float)gsl_vector_get(solver->x, 1);
  center[1] = (float)gsl_vector_get(solver->x, 2);
  *sigma    = 2.0 * (float)gsl_vector_get(solver->x, 3);
  *sigma    = 1.0 / (sqrt(*sigma)); /* back to sigma */

  /* cleanup */
  gsl_multifit_fdfsolver_free(solver);
  return data.RMS;
#else  /* No GSL - stubb */
  gchar *routine = "FitCGauss";
  Obit_log_error(err, OBIT_Error, 
		 "%s: GSL not available - cannot do fit", 
		     routine);
  return -1.0;
#endif /* GSL stuff */
} /* end FitCalGauss */

