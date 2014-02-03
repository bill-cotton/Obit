/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009-2014                                          */
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
/* Utility routine for fast sine/cosine calculation                   */
#include "ObitSinCos.h"
#include <math.h>

#define OBITSINCOSNTAB  1024  /* tabulated points per turn */
#define OBITSINCOSNTAB4 256   /* 1/4 size of table */
/** Is initialized? */
gboolean isInit = FALSE;
/** Sine lookup table covering 1 1/4 turn of phase */
ofloat sincostab[OBITSINCOSNTAB+OBITSINCOSNTAB4+1];
/** Angle spacing (radian) in table */
ofloat delta;
/** 1/2pi */
ofloat itwopi = 1.0/ (2.0 * G_PI);
/** 2pi */
ofloat twopi = (2.0 * G_PI);

/** AVX implementation 8 floats in parallel */
#if HAVE_AVX==1
#include <immintrin.h>

typedef __m256  v8sf;
typedef __m256i v4si;

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
  v4si      v;
} ALIGN32_END V4SI;

/* Constants */
#define _OBIT_TWOPI  (2.0 * G_PI)           /* 2pi */
#define _OBIT_ITWOPI 1.0/ (2.0 * G_PI)      /* 1/2pi */
#define _OBIT_DELTA  0.0061359231515425647  /* table spacing = 2pi/Obit_NTAB */
#define _OBIT_NTAB   OBITSINCOSNTAB         /* size of table -1 */
#define _OBIT_NTAB4  OBITSINCOSNTAB4
static const v8sf _half = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}; /* 0.5 vector */
static const v8sf _ntab = { _OBIT_NTAB, _OBIT_NTAB, _OBIT_NTAB, _OBIT_NTAB,
			    _OBIT_NTAB, _OBIT_NTAB, _OBIT_NTAB, _OBIT_NTAB};
static const v8sf _i2pi = { _OBIT_ITWOPI, _OBIT_ITWOPI, _OBIT_ITWOPI, _OBIT_ITWOPI, 
			    _OBIT_ITWOPI, _OBIT_ITWOPI, _OBIT_ITWOPI, _OBIT_ITWOPI};
static const v8sf _toopi= { _OBIT_TWOPI, _OBIT_TWOPI, _OBIT_TWOPI, _OBIT_TWOPI, 
			    _OBIT_TWOPI, _OBIT_TWOPI, _OBIT_TWOPI, _OBIT_TWOPI};
static const v8sf _dlta = { _OBIT_DELTA, _OBIT_DELTA, _OBIT_DELTA, _OBIT_DELTA, 
			    _OBIT_DELTA, _OBIT_DELTA, _OBIT_DELTA, _OBIT_DELTA};
/** 
 * Fast AVX (8) vector sine/cosine of angle 
 * Approximate sine/cosine, no range or value checking
 * \param angle  angle in radians
 * \param table  lookup table
 * \param s      [out] sine(angle)
 * \param c      [out] cosine(angle)
*/
void fast_sincos_ps(v8sf angle, float *table, v8sf *s, v8sf *c) {
  v8sf anglet, temp, ft, cell, it, sine, cosine, d;
  V4SI addr;

  /* get angle in turns */
  anglet = _mm256_mul_ps(angle, _i2pi);
  
  /* truncate to [0,1] turns */
  temp   = _mm256_floor_ps(anglet);        /* next lowest (signed) integer  */
  ft     = _mm256_sub_ps (anglet, temp);   /* Fractional turn */

  /* Table lookup, cos(phase) = sin(phase + 1/4 turn)*/
  it        = _mm256_mul_ps(ft, _ntab);    /* To cells in table */
  it        = _mm256_add_ps(it, _half);    /* add half */
  cell      = _mm256_floor_ps(it);         /* round to nearest cell */
  addr.v    = _mm256_cvtps_epi32(cell);    /* to integers */
  /* use union to load sine and cosine values */
  sine      = _mm256_set_ps(table[addr.i[7]],
			    table[addr.i[6]],
			    table[addr.i[5]],
			    table[addr.i[4]],
			    table[addr.i[3]],
			    table[addr.i[2]],
			    table[addr.i[1]],
			    table[addr.i[0]]);
  cosine    = _mm256_set_ps(table[addr.i[7]+_OBIT_NTAB4],
			    table[addr.i[6]+_OBIT_NTAB4],
			    table[addr.i[5]+_OBIT_NTAB4],
			    table[addr.i[4]+_OBIT_NTAB4],
			    table[addr.i[3]+_OBIT_NTAB4],
			    table[addr.i[2]+_OBIT_NTAB4],
			    table[addr.i[1]+_OBIT_NTAB4],
			    table[addr.i[0]+_OBIT_NTAB4]);

  /* One term Taylor series  */
  anglet = _mm256_mul_ps(ft, _toopi);       /* Now angle in radians [0,2 pi] */

  d      = _mm256_mul_ps (cell, _dlta);     /* tabulated phase = cell*delta_phase */
  d      = _mm256_sub_ps (anglet, d);       /* actual-tabulated phase */
  /* Cosine */
  temp   = _mm256_mul_ps (sine,d);
  *c     = _mm256_sub_ps (cosine, temp);
  /* Sine */
  temp   = _mm256_mul_ps (cosine,d);
  *s     = _mm256_add_ps (sine, temp);

  _mm_empty();  /* wait for operations to finish */
  return ;
} /* end AVX fast_sincos_ps */

/* end HAVE_AVX */

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

/* Constants */
#define _OBIT_TWOPI  6.2831853071795862     /* 2pi */
#define _OBIT_ITWOPI 0.15915494309189535    /* 1/2pi */
#define _OBIT_DELTA  0.0061359231515425647  /* table spacing = 2pi/Obit_NTAB */
#define _OBIT_NTAB   1024.0                 /* size of table -1 */
#define _OBIT_NTAB4  256
static const v4sf _half = {0.5, 0.5, 0.5, 0.5}; /* 0.5 vector */
static const v4sf _ntab = { _OBIT_NTAB, _OBIT_NTAB, _OBIT_NTAB, _OBIT_NTAB};
static const v4sf _i2pi = { _OBIT_ITWOPI, _OBIT_ITWOPI, _OBIT_ITWOPI, _OBIT_ITWOPI};
static const v4sf _toopi= { _OBIT_TWOPI, _OBIT_TWOPI, _OBIT_TWOPI, _OBIT_TWOPI};
static const v4sf _dlta = { _OBIT_DELTA, _OBIT_DELTA, _OBIT_DELTA, _OBIT_DELTA};

/** 
 * Fast SSE (4) vector sine/cosine of angle 
 * Approximate sine/cosine, no range or value checking
 * \param angle  angle in radians
 * \param table  lookup table
 * \param s      [out] sine(angle)
 * \param c      [out] cosine(angle)
*/
void fast_sincos_ps(v4sf angle, float *table, v4sf *s, v4sf *c) {
  v4sf anglet, temp, it, zero, mask, one, sine, cosine, d;
  v2si itLo, itHi;
  V2SI iaddrLo, iaddrHi;

  /* angle in turns */
  anglet = _mm_mul_ps(angle, _i2pi);
  
  /* truncate to [0,1] turns */
  /* Get full turns */
  itLo   = _mm_cvttps_pi32 (anglet);       /* first two truncated */
  temp   = _mm_movehl_ps (anglet,anglet);  /* upper two values into lower */
  itHi   = _mm_cvttps_pi32 (temp);         /* second two truncated */
  it     = _mm_cvtpi32_ps (temp, itHi);    /* float upper values */
  temp   = _mm_movelh_ps (it, it);         /* swap */
  it     = _mm_cvtpi32_ps (temp, itLo);    /* float lower values */

  /* If anglet negative, decrement it */
  zero   = _mm_setzero_ps ();              /* Zeros */
  mask   = _mm_cmplt_ps (anglet,zero);     /* Comparison to mask */
  one    = _mm_set_ps1 (1.0);              /* ones */
  one    = _mm_and_ps(one, mask);          /* mask out positive values */
  it     = _mm_sub_ps (it, one);
  anglet = _mm_sub_ps (anglet, it);        /* fold to [0,2pi] */

  /* Table lookup, cos(phase) = sin(phase + 1/4 turn)*/
  it        = _mm_mul_ps(anglet, _ntab);           /* To cells in table */
  it        = _mm_add_ps(it, _half);               /* To cells in table */
  iaddrLo.v = _mm_cvttps_pi32 (it);
  temp      = _mm_movehl_ps (it,it);               /* Round */
  iaddrHi.v = _mm_cvttps_pi32 (temp);
  sine      = _mm_setr_ps(table[iaddrLo.i[0]],table[iaddrLo.i[1]],
			  table[iaddrHi.i[0]],table[iaddrHi.i[1]]);
  cosine    = _mm_setr_ps(table[iaddrLo.i[0]+_OBIT_NTAB4],
			  table[iaddrLo.i[1]+_OBIT_NTAB4],
			  table[iaddrHi.i[0]+_OBIT_NTAB4],
			  table[iaddrHi.i[1]+_OBIT_NTAB4]);

  /* One term Taylor series  */
  anglet = _mm_mul_ps(anglet, _toopi);          /* Now angle in radians */
  temp   = _mm_cvtpi32_ps (it, iaddrHi.v);      /* float upper values */
  it     = _mm_movelh_ps (temp,temp);           /* swap */
  it     = _mm_cvtpi32_ps (it, iaddrLo.v);      /* float lower values */
  d      = _mm_mul_ps (it, _dlta);              /* tabulated phase */
  d      = _mm_sub_ps (anglet, d);              /* actual-tabulated phase */
  /* Cosine */
  temp   = _mm_mul_ps (sine,d);
  *c     = _mm_sub_ps (cosine, temp);
  /* Sine */
  temp   = _mm_mul_ps (cosine,d);
  *s     = _mm_add_ps (sine, temp);

  _mm_empty();  /* wait for operations to finish */
  return ;
} /* end SSE fast_sincos_ps */

#endif  /* HAVE_AVX */

/** 
 * Initialization 
 */
 void ObitSinCosInit(void)
{
  olong i;
  float angle;

  if (isInit) return;
  isInit = TRUE;  /* Now initialized */
  delta = ((2.0 *G_PI)/OBITSINCOSNTAB);

  for (i=0; i<(OBITSINCOSNTAB+OBITSINCOSNTAB4+1); i++) {
    angle = delta * i;
    sincostab[i] = (ofloat)sinf(angle);
  }
} /* end ObitSinCosInit */

/** 
 * Calculate sine/cosine of angle 
 * Lookup table initialized on first call
 * \param angle  angle in radians
 * \param sin    [out] sine(angle)
 * \param cos    [out] cosine(angle)
*/
void ObitSinCosCalc(ofloat angle, ofloat *sin, ofloat *cos)
{
  olong it, itt;
  ofloat anglet, ss, cc, d;

  /* Initialize? */
  if (!isInit) ObitSinCosInit();

  /* angle in turns */
  anglet = angle*itwopi;

  /* truncate to [0,1] turns */
  it = (olong)anglet;     
  if (anglet<0.0) it--;   /* fold to positive */
  anglet -= it;

  /* Lookup, cos(phase) = sin(phase + 1/4 turn) */
  itt = (olong)(0.5 + anglet*OBITSINCOSNTAB);
  ss  = sincostab[itt];
  cc  = sincostab[itt+OBITSINCOSNTAB4];

  /* One term Taylor series  */
  d = anglet*twopi - delta*itt;
  *sin = ss + cc * d;
  *cos = cc - ss * d;
} /* end ObitSinCosCalc */

/** 
 * Calculate sine/cosine of vector of angles, uses AVX or SSE implementation if available
 * Lookup table initialized on first call
 * \param n      Number of elements to process
 * \param angle  array of angles in radians
 * \param sin    [out] sine(angle)
 * \param cos    [out] cosine(angle)
*/
void ObitSinCosVec(olong n, ofloat *angle, ofloat *sin, ofloat *cos)
{
  olong i, nleft, it, itt;
  ofloat anglet, ss, cc, d;
  /** SSE implementation */
#if   HAVE_AVX==1
  olong ndo;
  v8sf vanglet, vss, vcc;
#elif HAVE_SSE==1
  olong ndo;
  V4SF vanglet, vss, vcc;
#endif /* HAVE_SSE */

  if (n<=0) return;
  
  /* Initialize? */
  if (!isInit) ObitSinCosInit();
  
  nleft = n;   /* Number left to do */
  i     = 0;   /* None done yet */

 /** avx implementation */
#if HAVE_AVX==1
  /* Loop in groups of 8 */
  ndo = nleft - nleft%8;  /* Only full groups of 8 */
  for (i=0; i<ndo; i+=8) {
    vanglet = _mm256_loadu_ps(angle); angle += 8;
      
    fast_sincos_ps(vanglet, sincostab, &vss, &vcc);
    _mm256_storeu_ps(sin, vss); sin += 8;
    _mm256_storeu_ps(cos, vcc); cos += 8;
 } /* end AVX loop */
 /** SSE implementation */
#elif HAVE_SSE==1
  /* Loop in groups of 4 */
  ndo = nleft - nleft%4;  /* Only full groups of 4 */
  for (i=0; i<ndo; i+=4) {
    vanglet.f[0] = *angle++;
    vanglet.f[1] = *angle++;
    vanglet.f[2] = *angle++;
    vanglet.f[3] = *angle++;    
    fast_sincos_ps(vanglet.v, sincostab, &vss.v, &vcc.v);
    *sin++ = vss.f[0];
    *sin++ = vss.f[1];
    *sin++ = vss.f[2];
    *sin++ = vss.f[3];
    *cos++ = vcc.f[0];
    *cos++ = vcc.f[1];
    *cos++ = vcc.f[2];
    *cos++ = vcc.f[3]; 
  } /* end SSE loop */
#endif /* HAVE_SSE */

  nleft = n-i;  /* How many left? */

 /* Loop doing any elements not done in AVX/SSE loop */
  for (i=0; i<nleft; i++) {
    /* angle in turns */
    anglet = (*angle++)*itwopi;
    
    /* truncate to [0,1] turns */
    it = (olong)anglet;     
    if (anglet<0.0) it--;   /* fold to positive */
    anglet -= it;
    
    /* Lookup, cos(phase) = sin(phase + 1/4 turn) */
    itt = (olong)(0.5 + anglet*OBITSINCOSNTAB);
    ss  = sincostab[itt];
    cc  = sincostab[itt+OBITSINCOSNTAB4];
    
    /* One term Taylor series  */
    d = anglet*twopi - delta*itt;
    *sin++ = ss + cc * d;
    *cos++ = cc - ss * d;
  } /* end loop over vector */
} /* end ObitSinCosVec */
