/* $Id$            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012                                               */
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
#include <glib.h>
#include "ObitTypes.h"
#ifndef OBITCOMPLEX_H 
#define OBITCOMPLEX_H 
/** 
 * \file ObitComplex.h
 * Obit complex arithmetic macro library.
 *
 * This allows complex arithmetic in C
 */
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/*--------------Class definitions-------------------------------------*/
/* Define complex structure */
typedef struct {
  ofloat real;
  ofloat imag;
} ocomplex;
typedef struct {
  odouble real;
  odouble imag;
} dcomplex;

/**
 * Set a complex value
 * \li out  = complex 
 * \li r    = real part
 * \li i    = imaginary part
 */
#define COMPLEX_SET(out, r, i) G_STMT_START{ \
  out.real = r; out.imag = i; \
}G_STMT_END  

/**
 * Set a complex exponential
 * out = exp(i arg)
 * \li out  = output 
 * \li arg  = argument in radians
 */
#define COMPLEX_EXP(out, arg)  G_STMT_START{ \
  out.real = cos(arg); out.imag = sin(arg);    \
}G_STMT_END  

/**
 * Add 2 complex values
 * out = in1 + in2
 * \li [out]out  = output complex
 * \li [in] in1  = input complex 
 * \li [in] in2  = input complex 
 */
#define COMPLEX_ADD2(out, in1, in2) G_STMT_START{ \
  out.real = in1.real + in2.real;   out.imag = in1.imag + in2.imag;  \
}G_STMT_END

/**
 * Add 3 complex values
 * out = in1 + in2 + in3
 * \li [out]out  = output complex
 * \li [in] in1  = input complex
 * \li [in] in2  = input complex
 * \li [in] in3  = input complex
 */
#define COMPLEX_ADD3(out, in1, in2, in3) G_STMT_START{ \
  out.real = in1.real + in2.real + in3.real;   \
  out.imag = in1.imag + in2.imag + in3.imag;  \
}G_STMT_END

/**
 * Add 4 complex values
 * out = in1 + in2 + in3 + in4
 * \li [out]out  = output complex
 * \li [in] in1  = input complex
 * \li [in] in2  = input complex
 * \li [in] in3  = input complex
 * \li [in] in4  = input complex
 */
#define COMPLEX_ADD4(out, in1, in2, in3, in4) G_STMT_START{ \
  out.real = in1.real + in2.real + in3.real + in4.real;   \
  out.imag = in1.imag + in2.imag + in3.imag + in4.imag;  \
}G_STMT_END

/**
 * Subtract complex values
 * out = in1 - in2
 * \li [out]out  = output complex
 * \li [in] in1  = input complex 
 * \li [in] in2  = input complex 
 */
#define COMPLEX_SUB(out, in1, in2)  G_STMT_START{ \
  out.real = in1.real - in2.real;   out.imag = in1.imag - in2.imag; \
}G_STMT_END

/**
 * Multiply 2 complex values
 * out = in1 * in2
 * \li [out]out  = output complex
 * \li [in] in1  = input complex 
 * \li [in] in2  = input complex 
 */
#define COMPLEX_MUL2(out, in1, in2)  G_STMT_START{      \
  out.real = in1.real * in2.real - in1.imag * in2.imag; \
  out.imag = in1.real * in2.imag + in1.imag * in2.real; \
}G_STMT_END

/**
 * Multiply 3 complex values
 * out = in1 * in2 * in3
 * \li [out]out  = output complex
 * \li [in] in1  = input complex 
 * \li [in] in2  = input complex 
 * \li [in] in3  = input complex 
 */
#define COMPLEX_MUL3(out, in1, in2, in3)  G_STMT_START{ \
{ ocomplex cmpxtmp1;                                    \
  cmpxtmp1.real = in1.real * in2.real - in1.imag * in2.imag; \
  cmpxtmp1.imag = in1.real * in2.imag + in1.imag * in2.real; \
  COMPLEX_MUL2 (out, cmpxtmp1, in3);                    \
}                                                       \
}G_STMT_END

/**
 * Multiply 4 complex values
 * out = in1 * in2 + in3 * in4
 * \li [out]out  = output complex
 * \li [in] in1  = input complex 
 * \li [in] in2  = input complex 
 * \li [in] in3  = input complex 
 * \li [in] in4  = input complex 
 */
#define COMPLEX_MUL4(out, in1, in2, in3, in4)  G_STMT_START{ \
{ ocomplex cmpxtmp1, cmpxtmp2;                               \
  cmpxtmp1.real = in1.real * in2.real - in1.imag * in2.imag; \
  cmpxtmp1.imag = in1.real * in2.imag + in1.imag * in2.real; \
  COMPLEX_MUL2 (cmpxtmp2, cmpxtmp1, in3);                    \
  COMPLEX_MUL2 (out, cmpxtmp2, in4);                         \
}                                                            \
}G_STMT_END

/**
 * Multiply 5 complex values
 * out = in1 * in2 * in3 * in4 * in5
 * \li [out]out  = output complex
 * \li [in] in1  = input complex 
 * \li [in] in2  = input complex 
 * \li [in] in3  = input complex 
 * \li [in] in4  = input complex 
 * \li [in] in5  = input complex 
 */
#define COMPLEX_MUL5(out, in1, in2, in3, in4, in5)  G_STMT_START{ \
{ ocomplex cmpxtmp1, cmpxtmp2;                                    \
  cmpxtmp1.real = in1.real * in2.real - in1.imag * in2.imag;      \
  cmpxtmp1.imag = in1.real * in2.imag + in1.imag * in2.real;      \
  COMPLEX_MUL2 (cmpxtmp2, cmpxtmp1, in3);                         \
  COMPLEX_MUL2 (cmpxtmp1, cmpxtmp2, in4);                         \
  COMPLEX_MUL2 (out, cmpxtmp1, in5);                              \
}                                                                 \
}G_STMT_END

/**
 * Multiply 6 complex values
 * out = in1 * in2 * in3 * in4 * in5 * in6
 * \li [out]out  = output complex
 * \li [in] in1  = input complex 
 * \li [in] in2  = input complex 
 * \li [in] in3  = input complex 
 * \li [in] in4  = input complex 
 * \li [in] in5  = input complex 
 * \li [in] in6  = input complex 
 */
#define COMPLEX_MUL6(out, in1, in2, in3, in4, in5, in6)  G_STMT_START{ \
{ ocomplex cmpxtmp1, cmpxtmp2;                                         \
  cmpxtmp1.real = in1.real * in2.real - in1.imag * in2.imag;           \
  cmpxtmp1.imag = in1.real * in2.imag + in1.imag * in2.real;           \
  COMPLEX_MUL2 (cmpxtmp2, cmpxtmp1, in3);                              \
  COMPLEX_MUL2 (cmpxtmp1, cmpxtmp2, in4);                              \
  COMPLEX_MUL2 (cmpxtmp2, cmpxtmp1, in5);                              \
  COMPLEX_MUL2 (out, cmpxtmp2, in6);                                   \
}                                                                      \
}G_STMT_END

/**
 * Multiply 7 complex values
 * out = in1 * in2 * in3 * in4 * in5 * in6 * in7
 * \li [out]out  = output complex
 * \li [in] in1  = input complex 
 * \li [in] in2  = input complex 
 * \li [in] in3  = input complex 
 * \li [in] in4  = input complex 
 * \li [in] in5  = input complex 
 * \li [in] in6  = input complex 
 * \li [in] in7  = input complex 
 */
#define COMPLEX_MUL7(out, in1, in2, in3, in4, in5, in6, in7)  G_STMT_START{ \
{ ocomplex cmpxtmp1, cmpxtmp2;                                         \
  cmpxtmp1.real = in1.real * in2.real - in1.imag * in2.imag;           \
  cmpxtmp1.imag = in1.real * in2.imag + in1.imag * in2.real;           \
  COMPLEX_MUL2 (cmpxtmp2, cmpxtmp1, in3);                              \
  COMPLEX_MUL2 (cmpxtmp1, cmpxtmp2, in4);                              \
  COMPLEX_MUL2 (cmpxtmp2, cmpxtmp1, in5);                              \
  COMPLEX_MUL2 (cmpxtmp1, cmpxtmp2, in6);                              \
  COMPLEX_MUL2 (out, cmpxtmp1, in7);                                   \
}                                                                      \
}G_STMT_END

/**
 * Multiply 8 complex values
 * out = in1 * in2 * in3 * in4 * in5 * in6 * in7 * in8
 * \li [out]out  = output complex
 * \li [in] in1  = input complex 
 * \li [in] in2  = input complex 
 * \li [in] in3  = input complex 
 * \li [in] in4  = input complex 
 * \li [in] in5  = input complex 
 * \li [in] in6  = input complex 
 * \li [in] in7  = input complex 
 * \li [in] in7  = input complex 
 */
#define COMPLEX_MUL8(out, in1, in2, in3, in4, in5, in6, in7, in8)  G_STMT_START{ \
{ ocomplex cmpxtmp1, cmpxtmp2;                                         \
  cmpxtmp1.real = in1.real * in2.real - in1.imag * in2.imag;           \
  cmpxtmp1.imag = in1.real * in2.imag + in1.imag * in2.real;           \
  COMPLEX_MUL2 (cmpxtmp2, cmpxtmp1, in3);                              \
  COMPLEX_MUL2 (cmpxtmp1, cmpxtmp2, in4);                              \
  COMPLEX_MUL2 (cmpxtmp2, cmpxtmp1, in5);                              \
  COMPLEX_MUL2 (cmpxtmp2, cmpxtmp2, in6);                              \
  COMPLEX_MUL2 (cmpxtmp2, cmpxtmp1, in7);                              \
  COMPLEX_MUL2 (out, cmpxtmp2, in8);                                   \
}                                                                      \
}G_STMT_END

/**
 * Conjugate a complex value
 * out = in *
 * \li [out]out  = output complex
 * \li [in] in   = input complex 
 */
#define COMPLEX_CONJUGATE(out, in) G_STMT_START{ \
  out.real = in.real;  out.imag = -in.imag; \
}G_STMT_END

/**
 * Negate a complex value
 * out = -in
 * \li [out]out  = output complex
 * \li [in] in   = input complex 
 */
#define COMPLEX_NEGATE(out, in) G_STMT_START{ \
  out.real = -in.real;  out.imag = -in.imag; \
}G_STMT_END

#endif /* OBITCOMPLEX_H */ 
