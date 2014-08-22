/* $Id: ObitExp.c 481 2014-05-28 19:40:41Z bill.cotton $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2011-2014                                          */
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
/* Utility routine for fast exp(-x) calculation                   */
#include "ObitVecFunc.h"
#include <math.h>

/** AVX implementation 8 floats in parallel */
#if HAVE_AVX==1
#include "avx_mathfun.h"
/** Natural log of array of 8 floats */
V8SF avx_log_ps(V8SF x) {
  return (V8SF) log256_ps((v8sf) x);
}
/** Exponential of array of 8 floats  */
V8SF avx_exp_ps(V8SF x) {
  return (V8SF) exp256_ps((v8sf) x);
}
/** Sine of array of 8 floats  */
V8SF avx_sin_ps(V8SF x) {
  return (V8SF) sin256_ps((v8sf) x);
}
/** Cosine of array of 8 floats  */
V8SF avx_cos_ps(V8SF x) {
  return (V8SF) cos256_ps((v8sf) x);
}
/** Sine and Cosine of array of 8 floats  */
void avx_sincos_ps(V8SF x, V8SF *s, V8SF *c) {
  sincos256_ps((v8sf) x, (v8sf*) s, (v8sf*) c);
}
/* end HAVE_AVX */

/** SSE implementation 4 floats in parallel */
#elif HAVE_SSE==1
#include  "sse_mathfun.h"
/** Natural log of array of 4 floats */
V4SF sse_log_ps(V4SF x) {
  return (V4SF) log_ps((v4sf) x)
}
/** Exponential of array of 4 floats  */
V4SF sse_exp_ps(V4SF x) {
  return (V4SF) exp_ps((v4sf) x);
}
/** Sine of array of 4 floats  */
V4SF sse_sin_ps(V4SF x) {
  return (V4SF) sin_ps((v4sf) x);
}
/** Cosine of array of 4 floats  */
V4SF sse_cos_ps(V4SF x) {
  return (V4SF) cos_ps((v4sf) x);
}
/** Sine and Cosine of array of 4 floats  */
void sse_sincos_ps(V4SF x, V4SF *s, V4SF *c) {
  sincos_ps((v4sf) x, (v4sf*) s, (v4sf*) c);
}
#endif  /* HAVE_SSE */
