/* $Id:  $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009                                               */
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

/** SSE implementation */
#ifdef HAVE_SSE
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
#define _PS_CONST(Name, Val)                                              \
  static const ALIGN16_BEG float _ps_##Name[4] ALIGN16_END = { Val, Val, Val, Val }
#define _PI32_CONST(Name, Val)                                            \
  static const ALIGN16_BEG int _pi32_##Name[4] ALIGN16_END = { Val, Val, Val, Val }

_PS_CONST(Obit_twopi,  6.2831853071795862);    /* 2pi */
_PS_CONST(Obit_itwopi,0.63661977236758138 );   /* 1/2pi */
_PS_CONST(Obit_delta, 0.0061359231515425647);  /* table spacing = 2pi/Obit_NTAB */
_PS_CONST(Obit_NTAB,  1024.0);                 /* size of table -1 */
_PS_CONST(Obit_HALF,  0.5);                    /* 0.5 */
_PI32_CONST(Obit_NTAB4,  256);                 /* 1/4 size of table */
#define _OBIT_TWOPI  6.2831853071795862     /* 2pi */
#define _OBIT_ITWOPI 0.15915494309189535    /* 1/2pi */
#define _OBIT_DELTA  0.0061359231515425647  /* table spacing = 2pi/Obit_NTAB */
#define _OBIT_NTAB   1024.0                 /* size of table -1 */
#define _OBIT_NTAB4  256

/** 
 * Fast sine/cosine of angle 
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
  temp   = _mm_set_ps1 (_OBIT_ITWOPI);
  anglet = _mm_mul_ps(angle, temp);
  
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
  temp      = _mm_set_ps1 (_OBIT_NTAB);
  it        = _mm_mul_ps(anglet, temp);           /* To cells in table */
  temp      = _mm_set_ps1 (0.5);
  it        = _mm_add_ps(it, temp);                   /* To cells in table */
  iaddrLo.v = _mm_cvttps_pi32 (it);
  temp      = _mm_movehl_ps (it,it);                  /* Round */
  iaddrHi.v = _mm_cvttps_pi32 (temp);
  sine      = _mm_setr_ps(table[iaddrLo.i[0]],table[iaddrLo.i[1]],
			  table[iaddrHi.i[0]],table[iaddrHi.i[1]]);
  cosine    = _mm_setr_ps(table[iaddrLo.i[0]+_OBIT_NTAB4],
			  table[iaddrLo.i[1]+_OBIT_NTAB4],
			  table[iaddrHi.i[0]+_OBIT_NTAB4],
			  table[iaddrHi.i[1]+_OBIT_NTAB4]);

  /* One term Taylor series  */
  temp   = _mm_set_ps1 (_OBIT_TWOPI);
  anglet = _mm_mul_ps(anglet, temp);            /* Now angle on radians */
  temp   = _mm_cvtpi32_ps (it, iaddrHi.v);      /* float upper values */
  it     = _mm_movelh_ps (temp,temp);           /* swap */
  it     = _mm_cvtpi32_ps (it, iaddrLo.v);      /* float lower values */
  temp   = _mm_set_ps1 (_OBIT_DELTA);
  d      = _mm_mul_ps (it, temp);               /* tabulated phase */
  d      = _mm_sub_ps (anglet, d);              /* actual-tabulated phase */
  /* Cosine */
  temp   = _mm_mul_ps (sine,d);
  *c     = _mm_sub_ps (cosine, temp);
  /* Sine */
  temp   = _mm_mul_ps (cosine,d);
  *s     = _mm_add_ps (sine, temp);

  _mm_empty();  /* wait for operations to finish */
  return ;
} /* end fast_sincos_ps */

#endif  /* HAVE_SSE */

/** 
 * Initialization 
 */
 void ObitSinCosInit(void)
{
  olong i;
  ofloat angle;

  isInit = TRUE;  /* Now initialized */
  delta = ((2.0 *G_PI)/OBITSINCOSNTAB);

  for (i=0; i<(OBITSINCOSNTAB+OBITSINCOSNTAB4+1); i++) {
    angle = delta * i;
    sincostab[i] = sinf(angle);
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
 * Calculate sine/cosine of vector of angles uses SSE implementation is available
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
#ifdef HAVE_SSE
  olong ndo;
  V4SF vanglet, vss, vcc;
#endif /* HAVE_SSE */
  
  /* Initialize? */
  if (!isInit) ObitSinCosInit();
  
  nleft = n;   /* Number left to do */
  i     = 0;   /* None done yet */

 /** SSE implementation */
#ifdef HAVE_SSE
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

 /* Loop doing any elements not done in SSE loop */
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
