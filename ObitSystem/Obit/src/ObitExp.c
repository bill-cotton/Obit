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
/* Utility routine for fast exp(-x) calculation                   */
#include "ObitExp.h"
#include <math.h>

#include "ObitVecFunc.h"
/** 
 * Init exp function - nop
 */
void ObitExpInit(void)
{
  return ; 

} /* end ObitExpInit */

/** 
 * Calculate exp of arg 
 * \param arg    argument
 * \return exp(arg)
 */
ofloat ObitExpCalc(ofloat arg)
{
  /* Fast scheme not accurate enough, use library */
  return expf(arg); 

} /* end ObitExpCalc */

/** 
 * Calculate exp(x) of vector of args uses AVX or SSE implementation if available
 * \param n      Number of elements to process
 * \param argarr array of args
 * \param exparr [out] exp(arg)
*/
void ObitExpVec(olong n, ofloat *argarr, ofloat *exparr)
{
  olong i, nleft;
#if   HAVE_AVX==1
  olong ndo;
  CV8SF varg, vexp;
#elif HAVE_SSE==1
  olong ndo;
  CV4SF vargt, vex;
#endif /* HAVE_SSE */
  
  if (n<=0) return;

  nleft = n;   /* Number left to do */
  i     = 0;   /* None done yet */

 /** avx implementation */
#if HAVE_AVX==1
  /* Loop in groups of 8 */
  ndo = nleft - nleft%8;  /* Only full groups of 8 */
  for (i=0; i<ndo; i+=8) {
    varg.v = _mm256_loadu_ps(argarr); argarr += 8;
    vexp.v = avx_exp_ps(varg.v);
    _mm256_storeu_ps(exparr, vexp.v); exparr += 8;
 } /* end AVX loop */
  /* Remainders, zero fill */
  nleft = n-i;  /* How many left? */
  for (i=0; i<nleft; i++) varg.f[i] = *argarr++;
  for (i=nleft; i<8; i++) varg.f[i] = 0.0;
  vexp.v = avx_exp_ps(varg.v);
  for (i=0; i<nleft; i++) *exparr++ = vexp.f[i];
 /** SSE implementation */
#elif HAVE_SSE==1
  /* Loop in groups of 4 */
  ndo = nleft - nleft%4;  /* Only full groups of 4 */
  for (i=0; i<ndo; i+=4) {
    vargt.f[0] = *argarr++;
    vargt.f[1] = *argarr++;
    vargt.f[2] = *argarr++;
    vargt.f[3] = *argarr++;
    vex.v = sse_exp_ps(vargt.v);
    *exparr++ = vex.f[0];
    *exparr++ = vex.f[1];
    *exparr++ = vex.f[2];
    *exparr++ = vex.f[3];
  } /* end SSE loop */
  /* Remainders, zero fill */
  nleft = n-i;  /* How many left? */
  for (i=0; i<nleft; i++) vargt.f[i] = *argarr++;
  for (i=nleft; i<4; i++) vargt.f[i] = 0.0;
  vex.v = sse_exp_ps(vargt.v);
  for (i=0; i<nleft; i++) *exparr++ = vex.f[i];
 /* end HAVE_SSE */
#else /* No SSE/AXV, use library function */
  for (i=0; i<n; i++) exparr[i] = expf(argarr[i]);
#endif
} /* end ObitExpVec */
