/* $Id: ObitFFT.c,v 1.7 2008/01/18 11:51:28 bcotton Exp $         */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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

#include <math.h>
#include "ObitFFT.h"
#include "ObitIOUVFITS.h"
#include "ObitIOUVAIPS.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFFT.c
 * ObitFFT class function definitions.
 * This class is derived from the Obit base class and is based on FFTW.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitFFT";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitFFTClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitFFTClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitFFTInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitFFTClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitFFTClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \param dir  The direction of the transform OBIT_FFT_Forward (R2C)
 *             or OBIT_FFT_Reverse (C2R).
 * \param type Whether OBIT_FFT_FullComplex (full C2C) or
 *             OBIT_FFT_HalfComplex (R2C or C2R).
 * \param rank of matrix range [1,7]
 * \param dim  dimensionality of each axis in column major (Fortran) order.
 *             If real/half complex is being used, then dim[0] should be
 *             the number of reals.
 * \return the new object.
 */
ObitFFT* newObitFFT (gchar* name, ObitFFTdir dir, ObitFFTtype type, 
		     olong rank, olong *dim)
{
  ObitFFT* out;
  olong i;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFFTClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitFFT));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");
  out->dir  = dir;
  out->type = type;
  out->rank = rank;
  g_assert (rank<=7); /* Not too many dimensions */
  for (i=0; i<rank; i++) out->dim[i] = dim[i];

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitFFTInit((gpointer)out);

 return out;
} /* end newObitFFT */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitFFTGetClass (void)
{
  return (gconstpointer)&myClassInfo;
} /* end ObitFFTGetClass */

/**
 * Given a length of data, suggest a larger size that will have an efficient FFT.
 * \param  length number of values to be transformed
 * \return a number equal or larger than length that will have an efficient 
 *  transform. 
 */
olong ObitFFTSuggestSize (gint length)
{
  olong i, out = length;
  /* Array of good values , even powers of 2, 3, 5 and small multiples */
  /* including odd ones  olong good[] = {
    2,3,4,5,6,8,9,10,12,16,
    20,24,25,27,32,40,48,64,80,81,
    96,125,128,160,192,243,256,300,320,384,400,432,480,512,576,600,
    625,640,768,800,864,960,1024,1280,1536,1600,1920,2048,2560,
    3072,3840,4096,5120,6144,6480,6561,8192,10240,12288,
    16384,19683,20480,24576,32768,40960,49152,65536
    81920,98304,131072,163840,177147,196608,234375,262144,
    327680,393216,524288,655360,786432,1048576,
    -1}; */
  /* really they should be even */
  olong good[] = {
    2,4,6,8,10,12,16,20,24,32,40,48,64,80,90,
    96,108,128,144,150,160,192,200,240,256,270,288,300,320,360,384,
    400,432,450,480,512,576,600,
    640,720,768,800,810,864,900,960,1024,1080,1152,
    1200,1280,1296,1350,1440,1500,1536,
    1600,1728,1800,1920,1944,2048,2160,2250,2400,2560,2592,2700,
    2880,3000,3072,3240,3456,3600,3840,4096,4320,4500,4800,
    5120,5400,5760,6000,6144,6480,6912,
    7500,7680,7776,8100,8192,8640,9000,
    10800,10240,12288,12960,13824,15552,16000,16384,
    17280,20480,20736,22500,24576,31104,32768,33750,40960,49152,65536,
    81920,98304,131072,163840,196608,262144,
    327680,393216,497664,524288,655360,786432,884736,
    1048576,1280000,1555200,1843200,2332800,2799360,
    -1};

  /* loop over acceptable values looking for next largest */
  i = 0;
  while (good[i]>0) {
    if (good[i]>=length) return good[i];
    i++;
  }

  /* No luck - return input value */
  return out;
} /* end ObitFFTSuggestSize */

/**
 * Do full real to half complex transform.
 * Must have been created with dir = OBIT_FFT_Forward and
 * type = OBIT_FFT_HalfComplex and have same geometry as constructor call.
 * \param in       Object with FFT structures.
 * \param inArray  Array to be transformed (undisturbed on output).
 * \param outArray Output array
 */
void ObitFFTR2C (ObitFFT *in, ObitFArray *inArray, ObitCArray *outArray)
{
#if HAVE_FFTW3==1   /* FFTW3 version */
#elif HAVE_FFTW==1  /* FFTW2 version */
#elif HAVE_GSL==1   /* GSL version */
  olong i, j, k, n, ni, nj, nk;
  size_t stride, dim;
  gboolean doFloat;
  gsl_complex_packed_array_float data;
  gsl_complex_packed_array dbldata;
  ofloat *idata=NULL, *odata=NULL;
#endif /* GSL */

  /* error checks */
  g_assert (ObitFFTIsA(in));
  g_assert (ObitFArrayIsA(inArray));
  g_assert (ObitCArrayIsA(outArray));
  g_assert (in->type == OBIT_FFT_HalfComplex);
  g_assert (in->dir  == OBIT_FFT_Forward);
  g_assert (in->rank == inArray->ndim);
  g_assert (in->dim[0] <= inArray->naxis[0]); /* Check two axes */
  g_assert (in->dim[1] <= inArray->naxis[1]); 
  g_assert (inArray->ndim == outArray->ndim);
  g_assert (inArray->naxis[0] <= 2*(outArray->naxis[0]-1));
  g_assert (inArray->naxis[1] <= outArray->naxis[1]);

  /* Lock ObitObjects aginst other threads */
  ObitThreadLock(in->thread);
  ObitThreadLock(inArray->thread);
  ObitThreadLock(outArray->thread);
   
  /* do transform */
#if HAVE_FFTW3==1   /* FFTW 3 */
  g_assert (in->RPlan!=NULL);
  fftwf_execute_dft_r2c (in->RPlan, (float*)inArray->array, 
		      (fftwf_complex*)outArray->array);

#elif HAVE_FFTW==1    /* FFTW 2 */
  g_assert (in->RPlan!=NULL);
  /* FFTW compiled in correct precision?*/
  g_assert (sizeof(fftw_real) == sizeof(ofloat));
  rfftwnd_one_real_to_complex (in->RPlan, (fftw_real*)inArray->array, 
			       (fftw_complex*)outArray->array);

#elif HAVE_GSL==1  /* Else try GSL version */

  /* Error checks */
  g_assert (in->RWavetab[0]!=NULL);
  if (in->rank>1) g_assert (in->CWavetab[1]!=NULL);
  if (in->rank>2) g_assert (in->CWavetab[2]!=NULL);
  g_assert (in->Rwork[0]!=NULL);
  if (in->rank>1) g_assert (in->Cwork[1]!=NULL);
  if (in->rank>2) g_assert (in->Cwork[2]!=NULL);

  /* Float or double? */
  doFloat = (sizeof(ofloat)==4);
  
  /* Loop over dimensions */
  ni = MAX(1, in->dim[1]);
  nj = MAX(1, 1+in->dim[0]/2);
  nk = ni*nj;
  /* First dimension */
  /* Use same convention as FFTW - half complex n/2+1 complex 
     FFT then unpack */
  idata = inArray->array;
  odata = outArray->array;
  stride = MAX (1, in->stride[0]);
  dim = MAX (1, in->dim[0]);
  for (i=0; i<ni; i++) {
    /* float or double? */
    if (doFloat) {
      data = (gsl_complex_packed_array_float)idata;
      gsl_fft_real_float_transform (data, stride, dim, in->RWavetab[0], in->Rwork[0]);
    } else { /* double */
      dbldata = (gsl_complex_packed_array)idata;
      gsl_fft_real_transform (dbldata, stride, dim, in->dblRWavetab[0], in->dblRwork[0]);
    }
    /* unpack to output */
    n = dim-1;
    *(odata)     = *(idata);
    *(odata+1)   = 0.0;
    *(odata+n+2) = 0.0;
    n *= sizeof(ofloat);
    memcpy (odata+2, idata+1, n);
    idata += dim;
    odata += dim+2;
  } /* end loop over first dimension */
  /* Second dimension if any */
  if (in->rank>1) {
    odata = outArray->array;
    stride = MAX (1, 1+in->stride[1]/2);
    dim = MAX (1, in->dim[1]);
    for (j=0; j<nj; j++) {
      /* float or double? */
      if (doFloat) {
	data = (gsl_complex_packed_array_float)odata;
	gsl_fft_complex_float_forward (data, stride, dim, in->CWavetab[1], in->Cwork[1]);
      } else { /* double */
	dbldata = (gsl_complex_packed_array)odata;
	gsl_fft_complex_forward (dbldata, stride, dim, in->dblCWavetab[1], in->dblCwork[1]);
      }
      odata += 2;
    } /* end loop over second dimension */
  } /* end rank 2 */
  /* Third dimension if any */
  if (in->rank>2) {
    odata = outArray->array;
    stride = MAX (1, in->stride[2]);
    dim = MAX (1, in->dim[2]);
    for (k=0; k<nk; k++) {
      /* float or double? */
      if (doFloat) {
	data = (gsl_complex_packed_array_float)odata;
	gsl_fft_complex_float_forward (data, stride, dim, in->CWavetab[2], in->Cwork[2]);
      } else { /* double */
	dbldata = (gsl_complex_packed_array)odata;
	gsl_fft_complex_forward (dbldata, stride, dim, in->dblCWavetab[2], in->dblCwork[0]);
      }
      odata += 2;
    } /* end loop over third dimension */
  } /* end rank 3 */
#endif /* FFT type */

  /* Unlock ObitObjects */
  ObitThreadUnlock(outArray->thread);
  ObitThreadUnlock(inArray->thread);
  ObitThreadUnlock(in->thread);
 
} /* end ObitFFTR2C */

/**
 * Do half complex to full real transform.
 * Must have been created with dir = OBIT_FFT_Reverse and
 * type = OBIT_FFT_HalfComplex and have same geometry as constructor call.
 * Note: FFT returned is not normalized.
 * \param in       Object with FFT structures.
 * \param inArray  Array to be transformed (disturbed on output).
 * \param outArray Output array
 */
void ObitFFTC2R (ObitFFT *in, ObitCArray *inArray, ObitFArray *outArray)
{
#if HAVE_FFTW3==1    /* FFTW 3 version */
#elif HAVE_FFTW==1    /* FFTW 2 version */
#elif HAVE_GSL==1   /* GSL version */
  olong i, j, k, n, ni, nj, nk;
  size_t stride, dim;
  gboolean doFloat;
  gsl_complex_packed_array_float data;
  gsl_complex_packed_array dbldata;
  ofloat *idata=NULL, *odata=NULL;
#endif /* GSL */

  /* error checks */
  g_assert (ObitFFTIsA(in));
  g_assert (ObitCArrayIsA(inArray));
  g_assert (ObitFArrayIsA(outArray));
  g_assert (in->type == OBIT_FFT_HalfComplex);
  g_assert (in->dir  == OBIT_FFT_Reverse);
  g_assert (in->rank == inArray->ndim);
  g_assert (in->dim[0] <= 2*(inArray->naxis[0]-1)); /* Check two axes */
  g_assert (in->dim[1] <= inArray->naxis[1]); 
  g_assert (inArray->ndim == outArray->ndim); 
  g_assert (2*(inArray->naxis[0]-1) == outArray->naxis[0]); 
  g_assert (inArray->naxis[1] == outArray->naxis[1]); 
 
  /* Lock ObitObjects aginst other threads */
  ObitThreadLock(in->thread);
  ObitThreadLock(inArray->thread);
  ObitThreadLock(outArray->thread);
   
  /* do transform */
#if HAVE_FFTW3==1   /* FFTW 3 */
  g_assert (in->RPlan!=NULL);
  fftwf_execute_dft_c2r (in->RPlan, (fftwf_complex*)inArray->array, 
     (float*)outArray->array);

#elif HAVE_FFTW==1     /* FFTW 2 version */
  g_assert (in->RPlan!=NULL);
  /* FFTW compiled in correct precision? */
  g_assert (sizeof(fftw_real) == sizeof(ofloat));
  rfftwnd_one_complex_to_real (in->RPlan, (fftw_complex*)inArray->array, 
			       (fftw_real*)outArray->array);
#elif HAVE_GSL==1  /* Else try GSL version */
  
  /* Error checks */
  g_assert (in->HCWavetab[0]!=NULL);
  if (in->rank>1) g_assert (in->CWavetab[1]!=NULL);
  if (in->rank>2) g_assert (in->CWavetab[2]!=NULL);
  g_assert (in->Rwork[0]!=NULL);
  if (in->rank>1) g_assert (in->Cwork[1]!=NULL);
  if (in->rank>2) g_assert (in->Cwork[2]!=NULL);
  
  /* Float or double? */
  doFloat = (sizeof(ofloat)==4);
  
  /* Loop over dimensions in reverse order */
  ni = MAX(1, in->dim[1]);
  nj = MAX(1, 1+in->dim[0]/2);
  nk = ni*nj;
  /* Third dimension if any */
  if (in->rank>2) {
    idata = inArray->array;
    stride = MAX (1, in->stride[2]);
    dim = MAX (1, in->dim[2]);
    for (k=0; k<nk; k++) {
      /* float or double? */
      if (doFloat) {
	data = (gsl_complex_packed_array_float)idata;
	gsl_fft_complex_float_backward (data, stride, dim, in->CWavetab[2], in->Cwork[2]);
      } else { /* double */
 	dbldata = (gsl_complex_packed_array)idata;
	gsl_fft_complex_backward (dbldata, stride, dim, in->dblCWavetab[2], in->dblCwork[2]);
      }
      idata += 2;
    } /* end loop over third dimension */
  } /* end rank 3 */
  /* Second dimension if any */
  if (in->rank>1) {
    idata = inArray->array;
    stride = MAX (1, 1+in->stride[1]/2);
    dim = MAX (1, in->dim[1]);
    for (j=0; j<nj; j++) {
      /* float or double? */
      if (doFloat) {
	data = (gsl_complex_packed_array_float)idata;
	gsl_fft_complex_float_backward (data, stride, dim, in->CWavetab[1], in->Cwork[1]);
      } else { /* double */
 	dbldata = (gsl_complex_packed_array)idata;
	gsl_fft_complex_backward (dbldata, stride, dim, in->dblCWavetab[1], in->dblCwork[1]);
      }
      idata += 2;
    } /* end loop over second dimension */
  } /* end rank 2 */
  /* First dimension */
  idata = inArray->array;
  odata = outArray->array;
  stride = MAX (1, in->stride[0]);
  dim = MAX (1, in->dim[0]);
  for (i=0; i<ni; i++) {
    /* Pack to output */
    n = (dim-1)*sizeof(ofloat);
    *(odata)     = *(idata);
    memcpy (odata+1, idata+2, n);
    /* float or double? */
    if (doFloat) {
      data = (gsl_complex_packed_array_float)odata;
      gsl_fft_halfcomplex_float_backward (data, stride, dim, in->HCWavetab[0], in->Rwork[0]);
    } else { /* double */
      dbldata = (gsl_complex_packed_array)odata;
      gsl_fft_halfcomplex_backward (dbldata, stride, dim, in->dblHCWavetab[0], in->dblRwork[0]);
    }
    idata += dim+2;
    odata += dim;
  } /* end loop over first dimension */
#endif /* HAVE_GSL */

  /* Unlock ObitObjects */
  ObitThreadUnlock(outArray->thread);
  ObitThreadUnlock(inArray->thread);
  ObitThreadUnlock(in->thread);
 } /* end ObitFFTC2R */

/**
 * Do full complex to complex transform in directions specified in in..
 * Must have been created with dir = OBIT_FFT_Reverse have same geometry as 
 * constructor call.
 * Transform is in the direction specified in constructor call.
 * \param in       Object with FFT structures.
 * \param inArray  Array to be transformed (disturbed on output).
 * \param outArray Output array
 */
void ObitFFTC2C (ObitFFT *in, ObitCArray *inArray, ObitCArray *outArray)
{
#if HAVE_FFTW3==1  /* FFTW version */
#elif HAVE_FFTW==1  /* FFTW version */
#elif HAVE_GSL==1   /* GSL version */
  olong i, j, k, ni, nj, nk;
  size_t stride, dim;
  gboolean doFloat;
  gsl_complex_packed_array_float data;
  gsl_complex_packed_array dbldata;
  ofloat *idata=NULL;
#endif /* GSL */

  /* error checks */
  g_assert (ObitFFTIsA(in));
  g_assert (ObitCArrayIsA(inArray));
  g_assert (ObitCArrayIsA(outArray));
  g_assert (in->type == OBIT_FFT_FullComplex);
  g_assert (in->rank == inArray->ndim);
  g_assert (in->dim[0] <= inArray->naxis[0]); /* Check two axes */
  g_assert (in->dim[1] <= inArray->naxis[1]); 
  g_assert (inArray->ndim == outArray->ndim);
  g_assert (inArray->naxis[0] <= outArray->naxis[0]);
  g_assert (inArray->naxis[1] <= outArray->naxis[1]);

  /* Lock ObitObjects aginst other threads */
  ObitThreadLock(in->thread);
  ObitThreadLock(inArray->thread);
  ObitThreadLock(outArray->thread);
   
  /* do transform */
#if HAVE_FFTW3==1   /* FFTW 3 */
  g_assert (in->CPlan!=NULL);
  fftwf_execute_dft (in->CPlan, (fftwf_complex*)inArray->array, 
		  (fftwf_complex*)outArray->array);

#elif HAVE_FFTW==1  /* FFTW 2 version */
  g_assert (in->CPlan!=NULL);
  /* FFTW compiled in correct precision? */
  g_assert(sizeof(fftw_real) == sizeof(ofloat));
  fftwnd_one (in->CPlan, (fftw_complex*)inArray->array, 
    (fftw_complex*)outArray->array);

#elif HAVE_GSL==1  /* Else try GSL version */

  /* Error check */
  for (i=0; i<in->rank; i++) {
    g_assert (in->CWavetab[i]!=NULL);
    g_assert (in->Cwork[i]!=NULL);
  }

  /* Float or double? */
  doFloat = (sizeof(ofloat)==4);
  
  /* First copy to output and work in place */
  ni = inArray->arraySize*sizeof(ofloat);
  memcpy (inArray->array, outArray->array, ni);
  /* Loop over dimensions */
  ni = MAX(1, in->dim[1]);
  nj = MAX(1, in->dim[0]);
  nk = ni*nj;

  /* Forward transform */
  if (in->dir==OBIT_FFT_Forward) {
    idata = outArray->array;
    stride = MAX (1, in->stride[0]);
    dim = MAX (1, in->dim[0]);
    for (i=0; i<ni; i++) {
      /* float or double? */
      if (doFloat) {
	data = (gsl_complex_packed_array_float)idata;
	gsl_fft_complex_float_forward (data, stride, dim, in->CWavetab[0], in->Cwork[0]);
      } else { /* double */
	dbldata = (gsl_complex_packed_array)idata;
	gsl_fft_complex_forward (dbldata, stride, dim, in->dblCWavetab[0], in->dblCwork[0]);
      }
      idata += 2*dim;
    } /* end loop over first dimension */
    if (in->rank>1) {
      idata = outArray->array;
      stride = MAX (1, in->stride[1]);
      dim = MAX (1, in->dim[1]);
      for (j=0; j<nj; j++) {
	/* float or double? */
	if (doFloat) {
	  data = (gsl_complex_packed_array_float)idata;
	  gsl_fft_complex_float_forward (data, stride, dim, in->CWavetab[1], in->Cwork[1]);
	} else { /* double */
	  dbldata = (gsl_complex_packed_array)idata;
	  gsl_fft_complex_forward (dbldata, stride, dim, in->dblCWavetab[1], in->dblCwork[1]);
	}
	idata += 2;
      } /* end loop over second dimension */
    } /* end rank 2 */
    if (in->rank>2) {
      idata = outArray->array;
      stride = MAX (1, in->stride[2]);
      dim = MAX (1, in->dim[2]);
      for (k=0; k<nk; k++) {
	/* float or double? */
	if (doFloat) {
	  data = (gsl_complex_packed_array_float)idata;
	  gsl_fft_complex_float_forward (data, stride, dim, in->CWavetab[2], in->Cwork[2]);
	} else { /* double */
	  dbldata = (gsl_complex_packed_array)idata;
	  gsl_fft_complex_forward (dbldata, stride, dim, in->dblCWavetab[2], in->dblCwork[2]);
	}
 	idata += 2;
     } /* end loop over third dimension */
    } /* end rank 3 */
    /* end forward */

  /* Reverse transform */
  } else if (in->dir==OBIT_FFT_Reverse) {
    idata = outArray->array;
    for (i=0; i<ni; i++) {
      stride = MAX (1, in->stride[0]);
      dim = MAX (1, in->dim[0]);
      /* float or double? */
      if (doFloat) {
	data = (gsl_complex_packed_array_float)idata;
	gsl_fft_complex_float_inverse (data, stride, dim, in->CWavetab[0], in->Cwork[0]);
      } else { /* double */
	dbldata = (gsl_complex_packed_array)idata;
	gsl_fft_complex_inverse (dbldata, stride, dim, in->dblCWavetab[0], in->dblCwork[0]);
      }
      idata += dim;
    } /* end loop over first dimension */
    if (in->rank>1) {
      idata = outArray->array;
      for (j=0; j<nj; j++) {
	stride = MAX (1, in->stride[1]);
	dim = MAX (1, in->dim[1]);
	/* float or double? */
	if (doFloat) {
	  data = (gsl_complex_packed_array_float)idata;
	  gsl_fft_complex_float_inverse (data, stride, dim, in->CWavetab[1], in->Cwork[1]);
	} else { /* double */
	  dbldata = (gsl_complex_packed_array)idata;
	  gsl_fft_complex_inverse (dbldata, stride, dim, in->dblCWavetab[1], in->dblCwork[1]);
	}
	idata += 2;
      } /* end loop over second dimension */
    } /* end rank 2 */
    if (in->rank>2) {
      idata = outArray->array;
      for (k=0; k<nk; k++) {
	stride = MAX (1, in->stride[2]);
	dim = MAX (1, in->dim[2]);
	/* float or double? */
	if (doFloat) {
	  data = (gsl_complex_packed_array_float)idata;
	  gsl_fft_complex_float_inverse (data, stride, dim, in->CWavetab[2], in->Cwork[2]);
	} else { /* double */
	  dbldata = (gsl_complex_packed_array)idata;
	  gsl_fft_complex_inverse (dbldata, stride, dim, in->dblCWavetab[2], in->dblCwork[2]);
	}
 	idata += 2;
     } /* end loop over third dimension */
    } /* end rank 3 */
  } /* end reverse */
#endif /* HAVE_GSL */

  /* Unlock ObitObjects */
  ObitThreadUnlock(outArray->thread);
  ObitThreadUnlock(inArray->thread);
  ObitThreadUnlock(in->thread);
 } /* end ObitFFTC2C */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitFFTClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitFFTClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitFFTClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitFFTClassInfoDefFn (gpointer inClass)
{
  ObitFFTClassInfo *theClass = (ObitFFTClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitFFTClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitFFTClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitFFTGetClass;
  theClass->newObit       = (newObitFP)newObitFFT;
  theClass->ObitCopy      = NULL;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitFFTClear;
  theClass->ObitInit      = (ObitInitFP)ObitFFTInit;
  /* New to this class */
  theClass->newObitFFT  = (newObitFFTFP)newObitFFT;
  theClass->ObitFFTR2C  = (ObitFFTR2CFP)ObitFFTR2C;
  theClass->ObitFFTC2R  = (ObitFFTC2RFP)ObitFFTC2R;
  theClass->ObitFFTC2C  = (ObitFFTC2CFP)ObitFFTC2C;

} /* end ObitFFTClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitFFTInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFFT *in = inn;
  olong i, dim[7];
#if HAVE_FFTW3==1   /* FFTW 3 version */
  unsigned flags;
  int fftwdim[10];
  int dir;
  float inArr[50], outArr[50];  /* Dummy arrays */
#elif HAVE_FFTW==1  /* FFTW 2 version */
  fftw_direction dir;
  olong flag;
  int fftwdim[10];

#elif HAVE_GSL==1  /* Else try GSL version */
  gboolean doFloat;
#endif /* HAVE_FFTW */

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread       = newObitThread();

#if HAVE_FFTW3==1 /* FFTW 3 version */
  in->CPlan        = NULL;
  in->RPlan        = NULL;

#elif HAVE_FFTW==1  /* FFTW 2 version */
  in->CPlan        = NULL;
  in->RPlan        = NULL;

#elif HAVE_GSL==1  /* Else try GSL version */
  for (i=0; i<3; i++) {
    in->CWavetab[i]  = NULL;
    in->HCWavetab[i] = NULL;
    in->RWavetab[i]  = NULL;
    in->Rwork[i]     = NULL;
    in->Cwork[i]     = NULL;
    in->dblCWavetab[i]  = NULL;
    in->dblHCWavetab[i] = NULL;
    in->dblRWavetab[i]  = NULL;
    in->dblRwork[i]     = NULL;
    in->dblCwork[i]     = NULL;
  }
#endif /* HAVE_FFTW */

  /* Reverse order of dimensions since FFTW uses row major and Obit
     uses column major */
  for (i=0; i<in->rank; i++) dim[in->rank-i-1] = in->dim[i];

  /* initialize FFT Plan */
#if HAVE_FFTW3==1  /* FFTW 3 version */
  /* DEBUG 
  fprintf (stderr, "DEBUG Using FFTW3 rank %d, dir %d, dim %d %d\n\n",
	   in->rank,in->dir, dim[0],dim[1]);*/

  /* Be sure ofloat set to float */
  g_assert (sizeof(ofloat)==sizeof(float));

  flags = FFTW_ESTIMATE; /* may only be done once - use cheap method */
  /* Note: NEVER set this to  FFTW_MEASURE as this will cause it to write over
     the dummy arrays inArr and outArr which it thinks are the size needed 
     for the FFT. */

  for (i=0; i<in->rank; i++) fftwdim[i] = (int)dim[i]; /* FFTW data type */
  if (in->type==OBIT_FFT_FullComplex) {
    if (in->dir==OBIT_FFT_Forward) dir = FFTW_FORWARD;
    else dir = FFTW_BACKWARD;
    if (in->rank==1) {
      in->CPlan = fftwf_plan_dft_1d(fftwdim[0],(fftwf_complex*)inArr, (fftwf_complex*)outArr,
				    dir, flags);
    } else if (in->rank==2) {
      in->CPlan = fftwf_plan_dft_2d(fftwdim[0], fftwdim[1], 
				    (fftwf_complex*)inArr, (fftwf_complex*)outArr,
				    dir, flags);
    } else {
      in->CPlan = fftwf_plan_dft((int)in->rank, fftwdim,
				 (fftwf_complex*)inArr,  (fftwf_complex*)outArr, 
				 dir, flags);
    }
  } else if (in->type==OBIT_FFT_HalfComplex) {
    /* Dimensions for FFTW3???
    fftwdim[in->rank-1] = fftwdim[in->rank-1]/2 + 1; */
    if (in->dir==OBIT_FFT_Forward) { /* R2C */
      if (in->rank==1) {
	in->RPlan = fftwf_plan_dft_r2c_1d(fftwdim[0],(float*)inArr, (fftwf_complex*)outArr,
					  flags);
      } else if (in->rank==2) {
	in->RPlan = fftwf_plan_dft_r2c_2d(fftwdim[0], fftwdim[1], 
					  (float*)inArr, (fftwf_complex*)outArr,
					  flags);
      } else {
	in->RPlan = fftwf_plan_dft_r2c((int)in->rank, fftwdim, 
				       (float*)inArr,  (fftwf_complex*)outArr, 
				       flags);
      }
    } else {  /* C2R */
      
      if (in->rank==1) {
	in->RPlan = fftwf_plan_dft_c2r_1d(fftwdim[0],(fftwf_complex*)inArr, (float*)outArr,
					  flags);
      } else if (in->rank==2) {
	in->RPlan = fftwf_plan_dft_c2r_2d(fftwdim[0], fftwdim[1], 
					  (fftwf_complex*)inArr, (float*)outArr,
					  flags);
      } else {
	in->RPlan = fftwf_plan_dft_c2r((int)in->rank, fftwdim,
				       (fftwf_complex*)inArr,  (float*)outArr, 
				       flags);
      }
    }
  }

#elif HAVE_FFTW==1  /* FFTW 2 version */
  flag = FFTW_ESTIMATE; /* may only be done once - use cheap method */
  for (i=0; i<in->rank; i++) fftwdim[i] = (int)dim[i]; /* FFTW data type */
  if (in->type==OBIT_FFT_FullComplex) {
    if (in->dir==OBIT_FFT_Forward) dir = FFTW_FORWARD;
    else dir = FFTW_BACKWARD;
    in->CPlan = fftwnd_create_plan((int)in->rank, fftwdim, dir, (int)flag);
  } else if (in->type==OBIT_FFT_HalfComplex) {
    if (in->dir==OBIT_FFT_Forward) dir = FFTW_REAL_TO_COMPLEX;
    else dir = FFTW_COMPLEX_TO_REAL;
    in->RPlan = rfftwnd_create_plan((int)in->rank, fftwdim, dir, (int)flag);
  }

#elif HAVE_GSL==1  /* Else try GSL version */
  g_assert (in->rank<=3);  /* Only up to rank 3 */

  /* Float or double? */
  doFloat = (sizeof(ofloat)==4);
  
  for (i=0; i<in->rank; i++) {
    /* Full complex or dimension >1 */
    if ((in->type==OBIT_FFT_FullComplex) || (i>0)) {
      /* Set strides */
      in->stride[i] = 1;
      if (i>=1) in->stride[i] = in->stride[i-1]* in->dim[i-1];
      /* Allocate work structures */
      if (doFloat) { /* float */
	in->CWavetab[i]  = gsl_fft_complex_wavetable_float_alloc(in->dim[i]);
	in->Cwork[i]     = gsl_fft_complex_workspace_float_alloc(in->dim[i]);
      } else { /* double */
	in->dblCWavetab[i]  = gsl_fft_complex_wavetable_alloc(in->dim[i]);
	in->dblCwork[i]     = gsl_fft_complex_workspace_alloc(in->dim[i]);
      }
    } else if (in->type==OBIT_FFT_HalfComplex) {
      /* Set strides */
      in->stride[i] = 1;
      if (i==1) in->stride[i] = in->stride[i-1]*(in->dim[i-1]/2+1);
      if (i>=2) in->stride[i] = in->stride[i-1]* in->dim[i-1];
      /* Allocate work structures */
      if (doFloat) { /* float */
	if (i==0) in->Rwork[i]     = gsl_fft_real_workspace_float_alloc(in->dim[i]);
	if (i>0)  in->Cwork[i]     = gsl_fft_complex_workspace_float_alloc(in->dim[i]);
      } else { /* double */
	if (i==0) in->dblRwork[i]     = gsl_fft_real_workspace_alloc(in->dim[i]);
	if (i>0)  in->dblCwork[i]     = gsl_fft_complex_workspace_alloc(in->dim[i]);
      }
      if (in->dir==OBIT_FFT_Forward) {
	/* Half complex forward (R to C) */
	if (doFloat) { /* float */
	  in->RWavetab[i]  = gsl_fft_real_wavetable_float_alloc(in->dim[i]);
	} else { /* double */
	  in->dblRWavetab[i]  = gsl_fft_real_wavetable_alloc(in->dim[i]);
	}
      } else {
	/* Half complex backwards (C to R) */
	if (doFloat) { /* float */
	  in->HCWavetab[i]  = gsl_fft_halfcomplex_wavetable_float_alloc(in->dim[i]);
	} else { /* double */
	  in->dblHCWavetab[i]  = gsl_fft_halfcomplex_wavetable_alloc(in->dim[i]);
	}
      }
    }
  }
#else  /* Die horrible death if attempt to use */
  g_error("ObitFFTInit: FFTW NOT implemented");
#endif /* HAVE_FFTW */
} /* end ObitFFTInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitFFT* cast to an Obit*.
 */
void ObitFFTClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFFT *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
#if HAVE_FFTW3==1  /* FFTW 3 version */
  i = 0;
  if (in->CPlan) fftwf_destroy_plan(in->CPlan); in->CPlan = NULL;
  if (in->RPlan) fftwf_destroy_plan(in->RPlan); in->RPlan = NULL;

#elif HAVE_FFTW==1  /* FFTW 2 version */
  i = 0;
  if (in->CPlan)  fftwnd_destroy_plan(in->CPlan); in->CPlan = NULL;
  if (in->RPlan) rfftwnd_destroy_plan(in->RPlan); in->RPlan = NULL;

#elif HAVE_GSL==1  /* Else try GSL version */
  for (i=0; i<in->rank; i++) {
    if (in->CWavetab[i]!=NULL)  gsl_fft_complex_wavetable_float_free(in->CWavetab[i]);     in->CWavetab[i]=NULL;
    if (in->HCWavetab[i]!=NULL) gsl_fft_halfcomplex_wavetable_float_free(in->HCWavetab[i]);in->HCWavetab[i]=NULL;
    if (in->RWavetab[i]!=NULL)  gsl_fft_real_wavetable_float_free(in->RWavetab[i]);        in->RWavetab[i]=NULL;
    if (in->Rwork[i]!=NULL)     gsl_fft_real_workspace_float_free (in->Rwork[i]);          in->Rwork[i]=NULL;
    if (in->Cwork[i]!=NULL)     gsl_fft_complex_workspace_float_free(in->Cwork[i]);        in->Cwork[i]=NULL;
    if (in->dblCWavetab[i]!=NULL)  gsl_fft_complex_wavetable_free(in->dblCWavetab[i]);     in->dblCWavetab[i]=NULL;
    if (in->dblHCWavetab[i]!=NULL) gsl_fft_halfcomplex_wavetable_free(in->dblHCWavetab[i]);in->dblHCWavetab[i]=NULL;
    if (in->dblRWavetab[i]!=NULL)  gsl_fft_real_wavetable_free(in->dblRWavetab[i]);        in->dblRWavetab[i]=NULL;
    if (in->dblRwork[i]!=NULL)     gsl_fft_real_workspace_free (in->dblRwork[i]);          in->dblRwork[i]=NULL;
    if (in->dblCwork[i]!=NULL)     gsl_fft_complex_workspace_free(in->dblCwork[i]);        in->dblCwork[i]=NULL;
  }
#endif /* HAVE_FFTW */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitFFTClear */

