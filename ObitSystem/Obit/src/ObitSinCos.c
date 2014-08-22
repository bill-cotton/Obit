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
#include "ObitVecFunc.h"
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
 * \param n      Number of elements to process
 * \param angle  array of angles in radians
 * \param sin    [out] sine(angle)
 * \param cos    [out] cosine(angle)
*/
void ObitSinCosVec(olong n, ofloat *angle, ofloat *sin, ofloat *cos)
{
  olong i, nleft;
  /** SSE implementation */
#if   HAVE_AVX==1
  olong ndo;
  CV8SF vanglet, vss, vcc;
#elif HAVE_SSE==1
  olong ndo;
  CV4SF vanglet, vss, vcc;
#endif /* HAVE_SSE */

  if (n<=0) return;
  nleft = n;
 /** AVX implementation */
#if HAVE_AVX==1
  /* Loop in groups of 8 */
  ndo = nleft - nleft%8;  /* Only full groups of 8 */
  for (i=0; i<ndo; i+=8) {
    vanglet.v = _mm256_loadu_ps(angle); angle += 8;     
    avx_sincos_ps(vanglet.v, &vss.v, &vcc.v);
    _mm256_storeu_ps(sin, vss.v); sin += 8;
    _mm256_storeu_ps(cos, vcc.v); cos += 8;
 } /* end AVX loop */
  /* Remainders, zero fill */
  nleft = n-i;  /* How many left? */
  for (i=0; i<nleft; i++) vanglet.f[i] = *angle++;
  for (i=nleft; i<8; i++) vanglet.f[i] = 0.0;
  avx_sincos_ps(vanglet.v, &vss.v, &vcc.v);
  for (i=0; i<nleft; i++) 
    {*sin++ = vss.f[i]; *cos++ = vcc.f[i]; }
 /** SSE implementation */
#elif HAVE_SSE==1
  /* Loop in groups of 4 */
  ndo = nleft - nleft%4;  /* Only full groups of 4 */
  for (i=0; i<ndo; i+=4) {
    vanglet.f[0] = *angle++;
    vanglet.f[1] = *angle++;
    vanglet.f[2] = *angle++;
    vanglet.f[3] = *angle++;    
    sse_sincos_ps(vanglet.v, &vss.v, &vcc.v);
    *sin++ = vss.f[0];
    *sin++ = vss.f[1];
    *sin++ = vss.f[2];
    *sin++ = vss.f[3];
    *cos++ = vcc.f[0];
    *cos++ = vcc.f[1];
    *cos++ = vcc.f[2];
    *cos++ = vcc.f[3]; 
  } /* end SSE loop */
  /* Remainders, zero fill */
  nleft = n-i;  /* How many left? */
  for (i=0; i<nleft; i++) vanglet.f[i] = *angle++;
  for (i=nleft; i<4; i++) vanglet.f[i] = 0.0;
  sse_sincos_ps(vanglet.v, &vss.v, &vcc.v);
  for (i=0; i<nleft; i++) 
    {*sin++ = vss.f[i]; *cos++ = vcc.f[i]; }
   /* end HAVE_SSE */
#else   /* Use table lookup if no SSE/AVX*/
  olong it, itt;
  ofloat anglet, ss, cc, d;
  /* Initialize? */
  if (!isInit) ObitSinCosInit();
  
  /* Loop using table lookup */
  for (i=0; i<n; i++) {
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
#endif
} /* end ObitSinCosVec */
