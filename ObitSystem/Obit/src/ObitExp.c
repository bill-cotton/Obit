/* $Id$ */
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
/* Utility routine for fast exp(x) calculation                        */
#include "ObitExp.h"
#include <math.h>

/* Coefficients of expansion for powers of x */
#define COEF0 1.0
#define COEF1 1.0/2.0
#define COEF2 1.0/(2.0*3.0)
#define COEF3 1.0/(2.0*3.0*4.0)
#define COEF4 1.0/(2.0*3.0*4.0*5.0)
#define COEF5 1.0/(2.0*3.0*4.0*5.0*6.0)
#define COEF6 1.0/(2.0*3.0*4.0*5.0*6.0*7.0)
#define COEF7 1.0/(2.0*3.0*4.0*5.0*6.0*7.0*8.0)

/** AVX implementation 8 floats in parallel */
#if HAVE_AVX==1
#include <immintrin.h>

typedef __m256  v8sf;
typedef __m256i v8si;

/* gcc or icc */
# define ALIGN32_BEG
# define ALIGN32_END __attribute__((aligned(32)))

/* Union allowing c interface */
typedef ALIGN32_BEG union {
  float f[8];
  int   i[8];
  v8sf   v;
} ALIGN32_END V8SF;

/* Union allowing c interface */
typedef ALIGN32_BEG union {
  int  i[8];
  v8si v;
} ALIGN32_END V8SI;

/* Constants */
static const v8sf _c0 = {COEF0,COEF0,COEF0,COEF0,COEF0,COEF0,COEF0,COEF0};
static const v8sf _c1 = {COEF1,COEF1,COEF1,COEF1,COEF1,COEF1,COEF1,COEF1};
static const v8sf _c2 = {COEF2,COEF2,COEF2,COEF2,COEF2,COEF2,COEF2,COEF2};
static const v8sf _c3 = {COEF3,COEF3,COEF3,COEF3,COEF3,COEF3,COEF3,COEF3};
static const v8sf _c4 = {COEF4,COEF4,COEF4,COEF4,COEF4,COEF4,COEF4,COEF4};
static const v8sf _c5 = {COEF5,COEF5,COEF5,COEF5,COEF5,COEF5,COEF5,COEF5};
static const v8sf _c6 = {COEF6,COEF6,COEF6,COEF6,COEF6,COEF6,COEF6,COEF6};
static const v8sf _c7 = {COEF7,COEF7,COEF7,COEF7,COEF7,COEF7,COEF7,COEF7};
/** 
 * Fast vector exp(arg) using AVX (8 float) instructions
 * \param arg    argument array
 * \param e      [out] array of exp(arg)
 */
void fast_exp_ps(v8sf arg,  v8sf *e) {
  v8sf sum, argp, tmp;

  /* Init sum, 0 and 1st power terms */
  sum   = _mm256_add_ps(_c0, arg);

  /* Init arg to power */
  argp  = _mm256_mul_ps(arg, arg);         /* arg^2 */
  tmp   = _mm256_mul_ps(argp, _c1);
  sum   = _mm256_add_ps(sum, tmp);

  argp  = _mm256_mul_ps(argp, arg);         /* arg^3 */
  tmp   = _mm256_mul_ps(argp, _c2);
  sum   = _mm256_add_ps(sum, tmp);

  argp  = _mm256_mul_ps(argp, arg);         /* arg^4 */
  tmp   = _mm256_mul_ps(argp, _c3);
  sum   = _mm256_add_ps(sum, tmp);

  argp  = _mm256_mul_ps(argp, arg);         /* arg^5 */
  tmp   = _mm256_mul_ps(argp, _c4);
  sum   = _mm256_add_ps(sum, tmp);

  argp  = _mm256_mul_ps(argp, arg);         /* arg^6 */
  tmp   = _mm256_mul_ps(argp, _c5);
  sum   = _mm256_add_ps(sum, tmp);

  argp  = _mm256_mul_ps(argp, arg);         /* arg^7 */
  tmp   = _mm256_mul_ps(argp, _c6);
  sum   = _mm256_add_ps(sum, tmp);

  argp  = _mm256_mul_ps(argp, arg);         /* arg^8 */
  tmp   = _mm256_mul_ps(argp, _c7);
  *e    = _mm256_add_ps(sum, tmp);

  return ;
} /* end fast_exp_ps */

/** SSE implementation 4 floats in parallel */
#elif HAVE_SSE==1
#include <xmmintrin.h>

typedef __m128 v4sf;
typedef __m64 v2si;

/* gcc or icc */
# define ALIGN16_BEG
# define ALIGN16_END __attribute__((aligned(16)))

/* Union allowing c interface */
typedef ALIGN16_BEG union {
  float f[4];
  int   i[4];
  v4sf  v;
} ALIGN16_END V4SF;

/* Union allowing c interface */
typedef ALIGN16_BEG union {
  int   i[2];
  v2si  v;
} ALIGN16_END V2SI;

/** 
 * Fast vector exp(arg) using SSE instructions
 * \param arg    argument array
 * \param e      [out] array of exp(arg)
 */
void fast_exp_ps(v4sf arg, v4sf *e) {
  v4sf sum, argp, tmp, coef;

  /* Init sum, 0 and 1st power terms */
  coef  = _mm_set_ps1 (COEF0);
  sum   = _mm_add_ps(coef, arg);

  /* Init arg to power */
  argp  = _mm_mul_ps(arg, arg);         /* arg^2 */
  coef  = _mm_set_ps1 (COEF1);
  tmp   = _mm_mul_ps(argp, coef);
  sum   = _mm_add_ps(sum, tmp);

  argp  = _mm_mul_ps(argp, arg);         /* arg^3 */
  coef  = _mm_set_ps1 (COEF2);
  tmp   = _mm_mul_ps(argp, coef);
  sum   = _mm_add_ps(sum, tmp);

  argp  = _mm_mul_ps(argp, arg);         /* arg^4 */
  coef  = _mm_set_ps1 (COEF3);
  tmp   = _mm_mul_ps(argp, coef);
  sum   = _mm_add_ps(sum, tmp);

  argp  = _mm_mul_ps(argp, arg);         /* arg^5 */
  coef  = _mm_set_ps1 (COEF4);
  tmp   = _mm_mul_ps(argp, coef);
  sum   = _mm_add_ps(sum, tmp);

  argp  = _mm_mul_ps(argp, arg);         /* arg^6 */
  coef  = _mm_set_ps1 (COEF5);
  tmp   = _mm_mul_ps(argp, coef);
  sum   = _mm_add_ps(sum, tmp);

  argp  = _mm_mul_ps(argp, arg);         /* arg^7 */
  coef  = _mm_set_ps1 (COEF6);
  tmp   = _mm_mul_ps(argp, coef);
  sum   = _mm_add_ps(sum, tmp);

  argp  = _mm_mul_ps(argp, arg);         /* arg^8 */
  coef  = _mm_set_ps1 (COEF7);
  tmp   = _mm_mul_ps(argp, coef);
  *e    = _mm_add_ps(sum, tmp);
  _mm_empty();  /* wait for operations to finish */
  return ;
} /* end fast_exp_ps */

#endif  /* HAVE_SSE */

/** 
 * Calculate exp of -arg 
 * arg<minTable => 1.0
 * arg>maxTable => 0.0
 * Lookup table initialized on first call
 * \param arg    argument
 * \return exp(arg)
*/
ofloat ObitExpCalc(ofloat arg)
{
  ofloat out, argp;
  return exp(arg);  /* library routine better optimized */
  out  = COEF0 + arg;
  argp = arg*arg;     /* arg^2 */
  out += argp * COEF1;
  argp = argp*arg;     /* arg^3 */
  out += argp * COEF2;
  argp = argp*arg;     /* arg^4 */
  out += argp * COEF3;
  argp = argp*arg;     /* arg^5 */
  out += argp * COEF4;
  argp = argp*arg;     /* arg^6 */
  out += argp * COEF5;
  argp = argp*arg;     /* arg^7 */
  out += argp * COEF6;
  argp = argp*arg;     /* arg^8 */
  out += argp * COEF7;
  return out;
} /* end ObitExpCalc */

/** 
 * Calculate exp(x) of vector of args uses SSE or AVX implementation if available
 * \param n      Number of elements to process
 * \param argarr array of args
 * \param exparr [out] exp(arg)
*/
void ObitExpVec(olong n, ofloat *argarr, ofloat *exparr)
{
  olong i, nleft;
  ofloat argt;
#if   HAVE_AVX==1
  olong ndo;
  v8sf varg, vexp;
#elif HAVE_SSE==1
  olong ndo;
  V4SF vargt, vex;
#endif /* HAVE_SSE */
  
  if (n<=0) return;

  nleft = n;   /* Number left to do */
  i     = 0;   /* None done yet */

 /** avx implementation */
#if HAVE_AVX==1
  /* Loop in groups of 8 */
  ndo = nleft - nleft%8;  /* Only full groups of 8 */
  for (i=0; i<ndo; i+=8) {
    varg = _mm256_loadu_ps(argarr); argarr += 8;
      
    fast_exp_ps(varg, &vexp);
    _mm256_storeu_ps(exparr, vexp); exparr += 8;
 } /* end AVX loop */
 /** SSE implementation */
#elif HAVE_SSE==1
  /* Loop in groups of 4 */
  ndo = nleft - nleft%4;  /* Only full groups of 4 */
  for (i=0; i<ndo; i+=4) {
    vargt.f[0] = *argarr++;
    vargt.f[1] = *argarr++;
    vargt.f[2] = *argarr++;
    vargt.f[3] = *argarr++;
    fast_exp_ps(vargt.v,  &vex.v);
    *exparr++ = vex.f[0];
    *exparr++ = vex.f[1];
    *exparr++ = vex.f[2];
    *exparr++ = vex.f[3];
  } /* end SSE loop */
#endif /* HAVE_SSE */

  nleft = n-i;  /* How many left? */

 /* Loop doing any elements not done in SSE/AVX loop */
  for (i=0; i<nleft; i++) {
    /* arg  */
    argt = (*argarr++);
    *exparr++ = ObitExpCalc(argt);
  } /* end loop over vector */
} /* end ObitExpVec */
