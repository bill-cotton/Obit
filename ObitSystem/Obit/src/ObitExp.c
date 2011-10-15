/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2011                                               */
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

#define OBITEXPTAB  512  /* tabulated points */
/** Is initialized? */
static gboolean isInit = FALSE;
/** Exp lookup table covering minTable to maxTable */
static ofloat exptab[OBITEXPTAB];
/** Arg spacing  in table */
static ofloat delta=0.02;
/** inverse of delta */
static ofloat idelta=50.0;
/** min value tabulated */
static ofloat minTable=1.0e-5;
/**max value tabulated */
static ofloat maxTable=10.0;

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

_PI32_CONST(Obit_NTAB,  512);     /* size of table */
_PS_CONST(Obit_delta, 0.02);      /* table spacing MUST match delta */
_PS_CONST(Obit_idelta, 50.0);     /* 1/table spacing */
_PS_CONST(Obit_minTable, 1.0e-5); /* minimum tabulated value, MUST match minTable */
_PS_CONST(Obit_maxTable, 10.0);   /* maximum tabulated value, MUST match maxTable */
_PS_CONST(Obit_HALF,  0.5);       /* 0.5 */
_PS_CONST(Obit_ONE,   1.0);       /* 1.0 */

#define _OBIT_DELTA     0.02   /* table spacing MUST match delta */
#define _OBIT_IDELTA    50.0   /* 1/table spacing  */
#define _OBIT_MINTABLE  1.0e-5 /* minimum tabulated value, MUST match minTable */
#define _OBIT_MAXTABLE  10.0   /* maximum tabulated value, MUST match maxTable */
#define _OBIT_HALF       0.5   /* 0.5 */
#define _OBIT_ONE        1.0   /* 1.0 */
#define _OBIT_NTAB       512   /* size of table */

/** 
 * Fast vector exp(-arg) using SSE instructions
 * \param arg    argument array
 * \param table  lookup table
 * \param e      [out] array of exp(-arg)
 */
void fast_exp_ps(v4sf arg, float *table, v4sf *e) {
  v4sf cellf, temp, it,  one, exptabl, d;
  V2SI iaddrLo, iaddrHi;

  /* Clip to range */
  temp  = _mm_set_ps1 (_OBIT_MINTABLE);
  arg   = _mm_max_ps (arg, temp);      /* Lower bound */
  temp  = _mm_set_ps1 (_OBIT_MAXTABLE); 
  arg   = _mm_min_ps (arg, temp);      /* Upper bound */

  /* get arg in cells */
  cellf = _mm_set_ps1 (_OBIT_MINTABLE);  /* Min table value */
  d     = _mm_sub_ps(arg, cellf);        /* arg-minTable */
  temp  = _mm_set_ps1 (_OBIT_IDELTA);    /* 1/table spacing */
  cellf = _mm_mul_ps(d, temp);           /* (arg-minTable)/table spacing */
  temp  = _mm_set_ps1 (_OBIT_HALF);      /* 0.5 */
  cellf = _mm_add_ps(cellf, temp);  
  iaddrLo.v = _mm_cvttps_pi32 (cellf);    /* Round lower half */
  temp      = _mm_movehl_ps (cellf,cellf);/* swap */
  iaddrHi.v = _mm_cvttps_pi32 (temp);     /* Round upper half */

  /* Fetch tabulated values */
  exptabl   = _mm_setr_ps(table[iaddrLo.i[0]],table[iaddrLo.i[1]],
			  table[iaddrHi.i[0]],table[iaddrHi.i[1]]);

  /* Get difference in arg from tabulated points */
  temp   = _mm_cvtpi32_ps (temp, iaddrHi.v);      /* float upper values */
  it     = _mm_movelh_ps (temp,temp);           /* swap */
  it     = _mm_cvtpi32_ps (it, iaddrLo.v);      /* float lower values */
  /* it now has the floated, truncated cells */
  temp  = _mm_set_ps1 (_OBIT_DELTA);            /* table spacing */
  temp  = _mm_mul_ps(it, temp);                 /* cell*delta */
  d     = _mm_sub_ps(d, temp);                  /* d = arg-cell*delta-minTable */

  /* One term Taylor's series */   
  one   = _mm_set_ps1 (1.0);                   /* ones */
  d     = _mm_sub_ps(one, d);                  /* 1-d */
  *e    = _mm_mul_ps(d, exptabl);              /* table[cell]*(1.0-d) */

  _mm_empty();  /* wait for operations to finish */
  return ;
} /* end fast_exp_ps */

#endif  /* HAVE_SSE */

/** 
 * Initialization 
 */
 void ObitExpInit(void)
{
  olong i, cell;
  ofloat arg;

  isInit = TRUE;  /* Now initialized */

  for (i=0; i<OBITEXPTAB-1; i++) {
    arg = minTable + delta * i;
    exptab[i] = exp(-arg);
  }
  exptab[OBITEXPTAB-1] = 0.0;  /* Higher values */

  /* Zero cell for maxTable */
  cell = (olong)((maxTable-minTable)*idelta + 0.5);
  exptab[cell] = 0.0;
} /* end ObitSinCosInit */

/** 
 * Calculate exp of -arg 
 * arg<minTable => 1.0
 * arg>maxTable => 0.0
 * Lookup table initialized on first call
 * \param arg    argument
 * \return exp(-arg)
*/
ofloat ObitExpCalc(ofloat arg)
{
  olong cell;
  ofloat out, d;

  /* Initialize? */
  if (!isInit) ObitExpInit();

  /* Range test */
  if (arg<minTable) return 1.0;
  if (arg>maxTable) return 0.0;


  /* Cell in lookup table */
  cell = (olong)((arg-minTable)*idelta + 0.5);

  /* Difference from tabulated value */
  d = arg - minTable -cell*delta;

  /* Lookup plus one Taylor term,  
     NB, d(exp(-x)/dx = -exp(-x) */
  out = exptab[cell]*(1.0-d);
  return out;
} /* end ObitExpCalc */

/** 
 * Calculate exp(-x) of vector of args uses SSE implementation is available
 * arg<minTable => 1.0
 * arg>maxTable => 0.0
 * Lookup table initialized on first call
 * \param n      Number of elements to process
 * \param argarr array of args
 * \param exparr [out] exp(-arg)
*/
void ObitExpVec(olong n, ofloat *argarr, ofloat *exparr)
{
  olong i, nleft, cell;
  ofloat argt, d;
  /** SSE implementation */
#ifdef HAVE_SSE
  olong ndo;
  V4SF vargt, vex;
#endif /* HAVE_SSE */
  
  /* Initialize? */
  if (!isInit) ObitExpInit();
  
  nleft = n;   /* Number left to do */
  i     = 0;   /* None done yet */

 /** SSE implementation */
#ifdef HAVE_SSE
  /* Loop in groups of 4 */
  ndo = nleft - nleft%4;  /* Only full groups of 4 */
  for (i=0; i<ndo; i+=4) {
    vargt.f[0] = *argarr++;
    vargt.f[1] = *argarr++;
    vargt.f[2] = *argarr++;
    vargt.f[3] = *argarr++;
    fast_exp_ps(vargt.v, exptab, &vex.v);
    *exparr++ = vex.f[0];
    *exparr++ = vex.f[1];
    *exparr++ = vex.f[2];
    *exparr++ = vex.f[3];
  } /* end SSE loop */
#endif /* HAVE_SSE */

  nleft = n-i;  /* How many left? */

 /* Loop doing any elements not done in SSE loop */
  for (i=0; i<nleft; i++) {
    /* arg  */
    argt = (*argarr++);

    /* range check */
    if (argt<minTable) {*exparr++=1.0; continue;}
    if (argt>maxTable) {*exparr++=0.0; continue;}

    /* Cell in lookup table */
    cell = (olong)((argt-minTable)*idelta + 0.5);
    
    /* Difference from tabulated value */
    d = argt - minTable -cell*delta;
    
    /* Lookup plus one Taylor term,  
       NB, d(exp(-x)/dx = -exp(-x) */
    *exparr++ = exptab[cell]*(1.0-d);
  } /* end loop over vector */
} /* end ObitExpVec */