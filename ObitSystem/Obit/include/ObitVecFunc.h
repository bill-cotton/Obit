/* $Id: $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014                                               */
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
/** 
 * Utility package for  access to 3rd party sse/avx vector functions
 * Avalibility depends on compiler settings HAVE_SSE and HAVE_AVX
 * as well as the actual implementation of sse and avx
 */

#ifndef OBITVECFUNC_H 
#define OBITVECFUNC_H 
/* gcc or icc */
# define ALIGN32_BEG
# define ALIGN32_END __attribute__((aligned(32)))

/** AVX implementation 8 floats in parallel */
#if HAVE_AVX==1
#include <immintrin.h>
/* Union allowing c interface */
typedef __m256  V8SF; // vector of 8 float (avx)
typedef ALIGN32_BEG union {
  float f[8];
  int   i[8];
  V8SF   v;
} ALIGN32_END CV8SF;
/** Natural log of array of 8 floats */
V8SF avx_log_ps(V8SF x);
/** Exponential of array of 8 floats  */
V8SF avx_exp_ps(V8SF x);
/** Sine of array of 8 floats  */
V8SF avx_sin_ps(V8SF x);
/** Cosine of array of 8 floats  */
V8SF avx_cos_ps(V8SF x);
/** Sine and Cosine of array of 8 floats  */
void avx_sincos_ps(V8SF x, V8SF *s, V8SF *c);
/* end HAVE_AVX */

/** SSE implementation 4 floats in parallel */
#elif HAVE_SSE==1
#include <xmmintrin.h>
typedef __m128 V4SF;  // vector of 4 float
typedef ALIGN16_BEG union {
  float f[4];
  int   i[4];
  V4SF  v;
} ALIGN16_END CV4SF;

/** Natural log of array of 4 floats */
V4SF sse_log_ps(V4SF x);
/** Exponential of array of 4 floats  */
V4SF sse_exp_ps(V4SF x);
/** Sine of array of 4 floats  */
V4SF sse_sin_ps(V4SF x);
/** Cosine of array of 4 floats  */
 V4SF sse_cos_ps(V4SF x);
/** Sine and Cosine of array of 4 floats  */
 void sse_sincos_ps(V4SF x, V4SF *s, V4SF *c);
#endif  /* HAVE_SSE */

#endif /* OBITVECFUNC_H */ 
