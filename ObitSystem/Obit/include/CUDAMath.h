/* $Id: $   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2026                                               */
/*;  Associated Universities, Inc. Washington DC, USA.                */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the CUDA ObitFArray structure       */
/*  This is intended to be included in a class structure definition   */
/** 
 * \file CUDAMath.h
 * Cuda complex arithmetic macro library.
 *
 * This allows complex arithmetic in C
 */
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
#include "ObitCUDAUtil.h"
#ifndef CUDAMATH_H 
#define CUDAMATH_H 

/*--------------Class definitions-------------------------------------*/
/* Define complex structure */
typedef struct {
  float real;
  float imag;
} cudaComplex;

/* Define Matx structure - only 2x2 complex supported */
/** data array - types set to array pointer*/
typedef struct {
  cudaComplex cpx[4];
} cudaMatx;

/*--------------- CUDA setup, GPU routines  ----------------*/
/* This is a CUDA routine */
#if HAVE_GPU==1  /* Have a GPU */
#if IS_CUDA==1   /* in cuda compiler? */
#include "ObitCUDAUtil.h"

/* OUTPUT *MUST* BE CALL BY REFERENCE */

/**
 * Set a complex value
 * \li out  = complex 
 * \li r    = real part
 * \li i    = imaginary part
 */
__device__ void CUDA_COMPLEX_SET(cudaComplex* out, float r, float i) { 
  out->real = r; out->imag = i; 
}   /* end CUDA_COMPLEX_SET*/

/**
 * Set a complex exponential
 * out = exp(i arg)
 * \li out  = output 
 * \li arg  = argument in radians
 */
__device__ void CUDA_COMPLEX_EXP(cudaComplex* out, float arg)  {
  __sincosf (arg, &out->imag, &out->real);
}  /* end CUDA_COMPLEX_EXP*/

/**
 * Get square of modulus
 * out = real*real+imag*imag
 * \li in   = complex value in question
 */
__device__ float CUDA_COMPLEX_MOD2(cudaComplex* in)  {
  float out;
  out = in->real*in->real+in->imag*in->imag;
  return out;
}  /*  end   CUDA_COMPLEX_MOD2*/

/**
 * Add 2 complex values
 * out = in1 + in2
 * \li [out]out  = output complex
 * \li [in] in1  = input complex 
 * \li [in] in2  = input complex 
 */
__device__ void CUDA_COMPLEX_ADD2(cudaComplex* out, cudaComplex in1,
				  cudaComplex in2) {
  out->real = in1.real + in2.real;   out->imag = in1.imag + in2.imag; 
} /* end CUDA_COMPLEX_ADD2 */

/**
 * Add 3 complex values
 * out = in1 + in2 + in3
 * \li [out]out  = output complex
 * \li [in] in1  = input complex
 * \li [in] in2  = input complex
 * \li [in] in3  = input complex
 */
__device__ void CUDA_COMPLEX_ADD3(cudaComplex* out, cudaComplex in1, 
				  cudaComplex in2,  cudaComplex in3) { 
  out->real = in1.real + in2.real + in3.real;  
  out->imag = in1.imag + in2.imag + in3.imag;
} /* end CUDA_COMPLEX_ADD3 */

/**
 * Add 4 complex values
 * out = in1 + in2 + in3 + in4
 * \li [out]out  = output complex
 * \li [in] in1  = input complex
 * \li [in] in2  = input complex
 * \li [in] in3  = input complex
 * \li [in] in4  = input complex
 */
__device__ void CUDA_COMPLEX_ADD4(cudaComplex* out, cudaComplex in1, 
				  cudaComplex in2, cudaComplex in3,
				  cudaComplex in4) { 
  out->real = in1.real + in2.real + in3.real + in4.real;
  out->imag = in1.imag + in2.imag + in3.imag + in4.imag; 
} /* end CUDA_COMPLEX_ADD4 */

/**
 * Subtract complex values
 * out = in1 - in2
 * \li [out]out  = output complex
 * \li [in] in1  = input complex 
 * \li [in] in2  = input complex 
 */
__device__ void CUDA_COMPLEX_SUB(cudaComplex* out, cudaComplex in1, cudaComplex in2)  { 
  out->real = in1.real - in2.real;   out->imag = in1.imag - in2.imag; 
} /* end CUDA_COMPLEX_SUB */

/**
 * Multiply 2 complex values
 * out = in1 * in2
 * \li [out]out  = output complex
 * \li [in] in1  = input complex 
 * \li [in] in2  = input complex 
 */
__device__ void CUDA_COMPLEX_MUL2(cudaComplex* out, cudaComplex in1, cudaComplex in2)  {
  float re = in1.real * in2.real - in1.imag * in2.imag;
  float im = in1.real * in2.imag + in1.imag * in2.real;
  out->real = re; out->imag = im; 
} /* end CUDA_COMPLEX_MUL2 */

/**
 * Multiply a complex value by the conjugate of another
 * out = in1 * conj(in2)
 * \li [out]out  = output complex
 * \li [in] in1  = 1st input complex 
 * \li [in] in2  = 2nd input complex 
 */
__device__ void CUDA_COMPLEX_MULC(cudaComplex* out, cudaComplex in1, cudaComplex in2)  {
  float re = +in1.real*in2.real + in1.imag*in2.imag;
  float im = -in1.real*in2.imag + in1.imag*in2.real;
  out->real = re; out->imag = im; 
} /* end CUDA_COMPLEX_MULC */

/**
 * Multiply 3 complex values
 * out = in1 * in2 * in3
 * \li [out]out  = output complex
 * \li [in] in1  = input complex 
 * \li [in] in2  = input complex 
 * \li [in] in3  = input complex 
 */
__device__ void CUDA_COMPLEX_MUL3(cudaComplex* out, cudaComplex in1,
				  cudaComplex in2,  cudaComplex in3)  {  
  cudaComplex cmpxtmp1;  
  cmpxtmp1.real = in1.real * in2.real - in1.imag * in2.imag;
  cmpxtmp1.imag = in1.real * in2.imag + in1.imag * in2.real;
  CUDA_COMPLEX_MUL2 (out, cmpxtmp1, in3);
} /* end CUDA_COMPLEX_MUL3 */

/**
 * Multiply 4 complex values
 * out = in1 * in2 + in3 * in4
 * \li [out]out  = output complex
 * \li [in] in1  = input complex 
 * \li [in] in2  = input complex 
 * \li [in] in3  = input complex 
 * \li [in] in4  = input complex 
 */
__device__ void CUDA_COMPLEX_MUL4(cudaComplex* out, cudaComplex in1,
				  cudaComplex in2,  cudaComplex in3,
				  cudaComplex in4)  {
  cudaComplex cmpxtmp1, cmpxtmp2;
  cmpxtmp1.real = in1.real * in2.real - in1.imag * in2.imag; 
  cmpxtmp1.imag = in1.real * in2.imag + in1.imag * in2.real; 
  CUDA_COMPLEX_MUL2 (&cmpxtmp2, cmpxtmp1, in3);
  CUDA_COMPLEX_MUL2 (out, cmpxtmp2, in4);
} /* end CUDA_COMPLEX_MUL4 */

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
__device__ void CUDA_COMPLEX_MUL5(cudaComplex* out, cudaComplex in1,
				  cudaComplex in2,  cudaComplex in3,
				  cudaComplex in4,  cudaComplex in5)  {
  cudaComplex cmpxtmp1, cmpxtmp2;
  cmpxtmp1.real = in1.real * in2.real - in1.imag * in2.imag;
  cmpxtmp1.imag = in1.real * in2.imag + in1.imag * in2.real;
  CUDA_COMPLEX_MUL2 (&cmpxtmp2, cmpxtmp1, in3);
  CUDA_COMPLEX_MUL2 (&cmpxtmp1, cmpxtmp2, in4);
  CUDA_COMPLEX_MUL2 (out, cmpxtmp1, in5);
} /* end CUDA_COMPLEX_MUL5 */

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
__device__ void CUDA_COMPLEX_MUL6(cudaComplex* out, cudaComplex in1,
				  cudaComplex in2,  cudaComplex in3,
				  cudaComplex in4,  cudaComplex in5,
				  cudaComplex in6)  {
  cudaComplex cmpxtmp1, cmpxtmp2;
  cmpxtmp1.real = in1.real * in2.real - in1.imag * in2.imag;
  cmpxtmp1.imag = in1.real * in2.imag + in1.imag * in2.real;
  CUDA_COMPLEX_MUL2 (&cmpxtmp2, cmpxtmp1, in3);
  CUDA_COMPLEX_MUL2 (&cmpxtmp1, cmpxtmp2, in4);
  CUDA_COMPLEX_MUL2 (&cmpxtmp2, cmpxtmp1, in5);
  CUDA_COMPLEX_MUL2 (out, cmpxtmp2, in6);
} /* end CUDA_COMPLEX_MUL6 */

/**
 * Divide complex values
 * out = in1 / in2
 * Returns 0 if mod(in2)=0
 * \li [out]out  = output complex
 * \li [in] in1  = input complex 
 * \li [in] in2  = input complex 
 */
__device__ void CUDA_COMPLEX_DIV(cudaComplex* out, cudaComplex in1,
				 cudaComplex in2)  { 
  double den = in2.real*in2.real + in2.imag*in2.imag;
  if (den>0.0) {
    den = 1.0 / den;
    out->real = den*(in1.real*in2.real + in1.imag*in2.imag);
    out->imag = den*(in1.imag*in2.real - in1.real*in2.imag);
  } else {
    out->real = 0.0; out->imag = 0.0;
  }
} /* end CUDA_COMPLEX_DIV */

/**
 * Conjugate a complex value
 * out = in *
 * \li [out]out  = output complex
 * \li [in] in   = input complex 
 */
__device__  void CUDA_COMPLEX_CONJUGATE(cudaComplex*out, cudaComplex in) { 
  out->real = in.real;  out->imag = -in.imag; 
}  /* end CUDA_COMPLEX_CONJUGATE */

/**
 * Negate a complex value
 * out = -in
 * \li [out]out  = output complex
 * \li [in] in   = input complex 
 */
__device__ void CUDA_COMPLEX_NEGATE(cudaComplex*out,cudaComplex in) { 
  out->real = -in.real;  out->imag = -in.imag; 
}  /* end CUDA_COMPLEX_NEGATE */

/**
 * Multiply a complex value by a scalar
 * out = in * scalar
 * \li [out]out  = output complex
 * \li [in] scalar= scalar to multiply
 * \li [in] in   = input complex 
 */
__device__ void CUDA_COMPLEX_SMULT(cudaComplex*out, float scalar, cudaComplex in) { 
  out->real = scalar*in.real;  out->imag = scalar*in.imag; 
}  /* end CUDA_COMPLEX_SMULT */

/************************ Matrices ***********************/
/* NO CAN DO */
/**
 * Test if a Matrix is blanked, checks 1st element
 * \li matx = cudaMatx structure modify.
 */
//__device__  int cudaMatxIsBlank(cudaMatx *matx) {
//  float fblank = CUDAMagicF();
//  int out = 0;
//  if (matx->cpx[0].real==fblank) out = 1;
//  return out;
//} /* end cudaMatxIsBlank */

/* Create 2x2 complex matrix, zero fill
Allocate/free from host 
__device__   cudaMatx cudaMatxCreate();
 MORE only pointer? */

/* delete a cudaMatx 
__device__   void cudaMatxFree (cudaMatx *in);
MORE only pointer?*/

/** Set 2x2  complex */
__device__ void cudaMatxZero2C(cudaMatx *matx) {
  matx->cpx[0].real = 0.0; matx->cpx[0].imag = 0.0; 
  matx->cpx[1].real = 0.0; matx->cpx[1].imag = 0.0; 
  matx->cpx[2].real = 0.0; matx->cpx[2].imag = 0.0; 
  matx->cpx[3].real = 0.0; matx->cpx[3].imag = 0.0; 
} /* end cudaMatxZero2C */

/** Set 2x2  complex */
__device__ void cudaMatxSet2C(cudaMatx *matx, 
			      float R00, float I00, float R01, float I01, 
			      float R10, float I10, float R11, float I11) {
  matx->cpx[0].real = R00; matx->cpx[0].imag = I00; 
  matx->cpx[1].real = R01; matx->cpx[1].imag = I01; 
  matx->cpx[2].real = R10; matx->cpx[2].imag = I10; 
  matx->cpx[3].real = R11; matx->cpx[3].imag = I11; 
} /* end cudaMatxSet2C */

/** Fetch 2x2  complex */
__device__ void cudaMatxGet2C(cudaMatx *matx, 
			      float *R00, float *I00, float *R01, float *I01, 
			      float *R10, float *I10, float *R11, float *I11) {
  *R00 = matx->cpx[0].real; *I00 = matx->cpx[0].imag;
  *R01 = matx->cpx[1].real; *I01 = matx->cpx[1].imag;
  *R10 = matx->cpx[2].real; *I10 = matx->cpx[2].imag;
  *R11 = matx->cpx[3].real; *I11 = matx->cpx[3].imag;
} /* end  cudaMatxGet2C */

/** Copy matrix */
__device__   cudaMatx* cudaMatxCopy(cudaMatx *in, cudaMatx *out) {
  /* Only 2x2 complex */
  out->cpx[0] = in->cpx[0];
  out->cpx[1] = in->cpx[1];
  out->cpx[2] = in->cpx[2];
  out->cpx[3] = in->cpx[3];
  return out;
} /* end cudaMatxCopy */

/** Zero values */
__device__ void cudaMatxZero(cudaMatx *out) {
  CUDA_COMPLEX_SET(&out->cpx[0], 0.0, 0.0);
  CUDA_COMPLEX_SET(&out->cpx[1], 0.0, 0.0);
  CUDA_COMPLEX_SET(&out->cpx[2], 0.0, 0.0);
  CUDA_COMPLEX_SET(&out->cpx[3], 0.0, 0.0);
} /* end cudaMatxZero */

/** Unit matrix */
__device__ void cudaMatxUnit(cudaMatx *out) {
  CUDA_COMPLEX_SET(&out->cpx[0], 1.0, 0.0);
  CUDA_COMPLEX_SET(&out->cpx[1], 0.0, 0.0);
  CUDA_COMPLEX_SET(&out->cpx[2], 0.0, 0.0);
  CUDA_COMPLEX_SET(&out->cpx[3], 1.0, 0.0);
} /* end cudaMatxUnit */

/** element by element Divide */
__device__   void cudaMatxDiv(cudaMatx *in1, cudaMatx *in2, cudaMatx *out) {
  CUDA_COMPLEX_DIV (&out->cpx[0], in1->cpx[0], in2->cpx[0]);
  CUDA_COMPLEX_DIV (&out->cpx[1], in1->cpx[1], in2->cpx[1]);
  CUDA_COMPLEX_DIV (&out->cpx[2], in1->cpx[2], in2->cpx[2]);
  CUDA_COMPLEX_DIV (&out->cpx[3], in1->cpx[3], in2->cpx[3]);
}  /* end cudaMatxDiv */

/** element by element Multiply */
__device__   void cudaMatxMult(cudaMatx *in1, cudaMatx *in2, cudaMatx *out) {
  cudaComplex tcpx;
 
  cudaMatxZero(out);  /* Zero output */

  CUDA_COMPLEX_MUL2(&tcpx, in1->cpx[0], in2->cpx[0]);
  CUDA_COMPLEX_ADD2(&out->cpx[0], out->cpx[0], tcpx);
  CUDA_COMPLEX_MUL2(&tcpx, in1->cpx[1], in2->cpx[2]);
  CUDA_COMPLEX_ADD2(&out->cpx[0], out->cpx[0], tcpx);

  CUDA_COMPLEX_MUL2(&tcpx, in1->cpx[0], in2->cpx[1]);
  CUDA_COMPLEX_ADD2(&out->cpx[1], out->cpx[1], tcpx);
  CUDA_COMPLEX_MUL2(&tcpx, in1->cpx[1], in2->cpx[3]);
  CUDA_COMPLEX_ADD2(&out->cpx[1], out->cpx[1], tcpx);

  CUDA_COMPLEX_MUL2(&tcpx, in1->cpx[2], in2->cpx[0]);
  CUDA_COMPLEX_ADD2(&out->cpx[2], out->cpx[2], tcpx);
  CUDA_COMPLEX_MUL2(&tcpx, in1->cpx[3], in2->cpx[2]);
  CUDA_COMPLEX_ADD2(&out->cpx[2], out->cpx[2], tcpx);

  CUDA_COMPLEX_MUL2(&tcpx, in1->cpx[2], in2->cpx[1]);
  CUDA_COMPLEX_ADD2(&out->cpx[3], out->cpx[3], tcpx);
  CUDA_COMPLEX_MUL2(&tcpx, in1->cpx[3], in2->cpx[3]);
  CUDA_COMPLEX_ADD2(&out->cpx[3], out->cpx[3], tcpx);
}  /* end cudaMatxMult */

/** Multiply by conjugate transpose */
__device__   void cudaMatxMultCT(cudaMatx *in1, cudaMatx *in2, cudaMatx *out) {
  cudaComplex tcpx, tcpx2;
 
  cudaMatxZero(out);  /* Zero output */

  CUDA_COMPLEX_CONJUGATE(&tcpx2, in2->cpx[0]);
  CUDA_COMPLEX_MUL2(&tcpx, in1->cpx[0], tcpx2);
  CUDA_COMPLEX_ADD2(&out->cpx[0], out->cpx[0], tcpx);
  CUDA_COMPLEX_CONJUGATE(&tcpx2, in2->cpx[1]);
  CUDA_COMPLEX_MUL2(&tcpx, in1->cpx[1], tcpx2);
  CUDA_COMPLEX_ADD2(&out->cpx[0], out->cpx[0], tcpx);

  CUDA_COMPLEX_CONJUGATE(&tcpx2, in2->cpx[2]);
  CUDA_COMPLEX_MUL2(&tcpx, in1->cpx[0], tcpx2);
  CUDA_COMPLEX_ADD2(&out->cpx[1], out->cpx[0], tcpx);
  CUDA_COMPLEX_CONJUGATE(&tcpx2, in2->cpx[3]);
  CUDA_COMPLEX_MUL2(&tcpx, in1->cpx[1], tcpx2);
  CUDA_COMPLEX_ADD2(&out->cpx[1], out->cpx[0], tcpx);

  CUDA_COMPLEX_CONJUGATE(&tcpx2, in2->cpx[0]);
  CUDA_COMPLEX_MUL2(&tcpx, in1->cpx[2], tcpx2);
  CUDA_COMPLEX_ADD2(&out->cpx[2], out->cpx[0], tcpx);
  CUDA_COMPLEX_CONJUGATE(&tcpx2, in2->cpx[1]);
  CUDA_COMPLEX_MUL2(&tcpx, in1->cpx[3], tcpx2);
  CUDA_COMPLEX_ADD2(&out->cpx[2], out->cpx[0], tcpx);

  CUDA_COMPLEX_CONJUGATE(&tcpx2, in2->cpx[2]);
  CUDA_COMPLEX_MUL2(&tcpx, in1->cpx[2], tcpx2);
  CUDA_COMPLEX_ADD2(&out->cpx[3], out->cpx[0], tcpx);
  CUDA_COMPLEX_CONJUGATE(&tcpx2, in2->cpx[3]);
  CUDA_COMPLEX_MUL2(&tcpx, in1->cpx[3], tcpx2);
  CUDA_COMPLEX_ADD2(&out->cpx[3], out->cpx[0], tcpx);
} /* end cudaMatxMultCT */

/** Add */
__device__ void cudaMatxAdd(cudaMatx *in1, cudaMatx *in2, cudaMatx *out) {
  /* Only 2x2 complex */
  CUDA_COMPLEX_ADD2(&out->cpx[0], in1->cpx[0], in2->cpx[0]);
  CUDA_COMPLEX_ADD2(&out->cpx[1], in1->cpx[1], in2->cpx[1]);
  CUDA_COMPLEX_ADD2(&out->cpx[2], in1->cpx[2], in2->cpx[2]);
  CUDA_COMPLEX_ADD2(&out->cpx[3], in1->cpx[3], in2->cpx[3]);
} /* end cudaMatxAdd */

/** Subtract in1-in2 */
__device__ void cudaMatxSub(cudaMatx *in1, cudaMatx *in2, cudaMatx *out) {
  /* Only 2x2 complex */
  CUDA_COMPLEX_SUB(&out->cpx[0], in1->cpx[0], in2->cpx[0]);
  CUDA_COMPLEX_SUB(&out->cpx[1], in1->cpx[1], in2->cpx[1]);
  CUDA_COMPLEX_SUB(&out->cpx[2], in1->cpx[2], in2->cpx[2]);
  CUDA_COMPLEX_SUB(&out->cpx[3], in1->cpx[3], in2->cpx[3]);
} /* end  cudaMatxSub */

/** Conjugate transpose */
__device__ void cudaMatxCTrans(cudaMatx *in, cudaMatx *out) {
  out->cpx[0] = in->cpx[0]; out->cpx[0].imag = -out->cpx[0].imag;
  out->cpx[1] = in->cpx[2]; out->cpx[1].imag = -out->cpx[1].imag;
  out->cpx[2] = in->cpx[1]; out->cpx[2].imag = -out->cpx[2].imag;
  out->cpx[3] = in->cpx[3]; out->cpx[3].imag = -out->cpx[3].imag;
} /* end CTrans */

/** Inverse perfect linear feed Jones matrix */
__device__ void cudaMatxIPerfLinJones(cudaMatx *out) {
  float elp_x=0.0, elp_y=0.0, ori_x=0.0, ori_y=1.57079632679;
  float angle[4], sina[4], cosa[4], Jones[8] ,Det[2], d;
  /*  Check that 2x2 complex? */

  angle[0] = 0.78539816+elp_x; angle[1] = 0.78539816-elp_y;
  angle[2] = ori_x;            angle[3] = ori_y;
  __sincosf (angle[0], &sina[0], &cosa[0]);
  __sincosf (angle[1], &sina[1], &cosa[1]);
  __sincosf (angle[2], &sina[2], &cosa[2]);
  __sincosf (angle[3], &sina[3], &cosa[3]);
  Jones[0] =  cosa[0]*cosa[2]; Jones[1] = -cosa[0]*sina[2];
  Jones[2] =  sina[0]*cosa[2]; Jones[3] =  sina[0]*sina[2];
  Jones[4] =  sina[1]*cosa[3]; Jones[5] = -sina[1]*sina[3];
  Jones[6] =  cosa[1]*cosa[3]; Jones[7] =  cosa[1]*sina[3];
  /* inverse of determinant */
  Det[0] = (Jones[0]*Jones[6] - Jones[1]*Jones[7]) - (Jones[2]*Jones[4] - Jones[3]*Jones[5]);
  Det[1] = (Jones[0]*Jones[7] + Jones[1]*Jones[6]) - (Jones[2]*Jones[5] + Jones[3]*Jones[4]);
  /* Inverse of determinant */
  d = Det[0]*Det[0] + Det[1]*Det[1];
  if (d!=0.0) d = 1.0 / d;
  else d = 1.0;
  Det[0] *=  d;
  Det[1] *= -d;
  /* Inverse matrix out */
  CUDA_COMPLEX_SET (&out->cpx[3],   Jones[0]*Det[0]-Jones[1]*Det[1],    Jones[0]*Det[1]+Jones[1]*Det[0]);
  CUDA_COMPLEX_SET (&out->cpx[2], -(Jones[4]*Det[0]-Jones[5]*Det[1]), -(Jones[4]*Det[1]+Jones[5]*Det[0]));
  CUDA_COMPLEX_SET (&out->cpx[1], -(Jones[2]*Det[0]-Jones[3]*Det[1]), -(Jones[2]*Det[1]+Jones[3]*Det[0]));
  CUDA_COMPLEX_SET (&out->cpx[0],   Jones[6]*Det[0]-Jones[7]*Det[1],    Jones[6]*Det[1]+Jones[7]*Det[0]);
} /* end cudaMatxIPerfLinJones */

/** Inverse perfect circular feed Jones matrix 
    Really just the perfect linear Jones matrix */
__device__ void cudaMatxIPerfCirJones(cudaMatx *out) {
  float elp_x=0.0, elp_y=0.0, ori_x=0.0, ori_y=1.57079632679;
  float angle[4], sina[4], cosa[4], Jones[8];

  angle[0] = 0.78539816+elp_x; angle[1] =  0.78539816-elp_y;
  angle[2] = ori_x;            angle[3] = ori_y;
  __sincosf (angle[0], &sina[0], &cosa[0]);
  __sincosf (angle[1], &sina[1], &cosa[1]);
  __sincosf (angle[2], &sina[2], &cosa[2]);
  __sincosf (angle[3], &sina[3], &cosa[3]);
  Jones[0] =  cosa[0]*cosa[2]; Jones[1] = -cosa[0]*sina[2];
  Jones[2] =  sina[0]*cosa[2]; Jones[3] =  sina[0]*sina[2];
  Jones[4] =  sina[1]*cosa[3]; Jones[5] = -sina[1]*sina[3];
  Jones[6] =  cosa[1]*cosa[3]; Jones[7] =  cosa[1]*sina[3];
  
  /* Matrix out */
  CUDA_COMPLEX_SET (&out->cpx[0], Jones[0], Jones[1]);
  CUDA_COMPLEX_SET (&out->cpx[1], Jones[2], Jones[3]);
  CUDA_COMPLEX_SET (&out->cpx[2], Jones[4], Jones[5]);
  CUDA_COMPLEX_SET (&out->cpx[3], Jones[6], Jones[7]); 
} /* end  cudaMatxIPerfCirJones */

/* Multiply a scalar times an cudaMatx */
__device__   void cudaMatxSMult (cudaMatx* in, float scalar, cudaMatx* out) {
  out->cpx[0].real = scalar * in->cpx[0].real; out->cpx[0].imag = scalar * in->cpx[0].imag; 
  out->cpx[1].real = scalar * in->cpx[1].real; out->cpx[1].imag = scalar * in->cpx[1].imag; 
  out->cpx[2].real = scalar * in->cpx[2].real; out->cpx[2].imag = scalar * in->cpx[2].imag; 
  out->cpx[3].real = scalar * in->cpx[3].real; out->cpx[3].imag = scalar * in->cpx[3].imag; 
} /* end cudaMatxSMult */

#endif // GPU code
#endif // CUDA code
#endif // CUDAMATH_H
