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
/*  Functions to correct visibilities for antenna pattern effects     */
/** 
 * \file CUDALinCorSum,h
 * Cuda fn to calculate model visibility including beam corrections
 * Linear feeds
 *
 */
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
#ifndef CUDALINCORSUM_H 
#define CUDALINCORSUM_H 
//damn #if IS_CUDA==1
#include "CUDAMath.h"
/*--------------- CUDA setup, GPU kernal  ----------------*/
/* HAVE to include source here to get it in the same fule as kernal */
/* This is a CUDA routine */
/** Public: Jones correct and sum SkyModel components, Lin pol feed, Stokes I 
 * Converts SkyModel comps to pseudovisibilities, multiply by Jones matrices and sum.
 * \param numComp  Number of SkyModel components
 * \param compFlux Component amplitude with any factors multiplied
 * \param sinArr   Array of sine of phase to apply, per numComp
 * \param cosArr   Array of cosine of phase to apply, per numComp
 * \param Jones1   Jones matrix array of first antenna of baseline, per numComp
 * \param Jones2   Jones matrix array of second antenna of baseline, per numComp
 * \param c2p      cosine (Parallactic angle*2)
 * \param s2p      sine (Parallactic angle*2)
 * \param sumVis   [out] 2x2 complex matrix for resultant sum
 */
 extern "C"
__device__ void cudaJonesCorSum_Lin_I (long numComp, float *compFlux, 
	                               float *sinArr, float *cosArr,
		                       cudaMatx **Jones1, cudaMatx **Jones2, 
				       float c2p, float s2p,
				       cudaMatx *sumVis)
{
  long i;
  cudaComplex v, t1, t2, t3, txx, tyy, txy, tyx;
  cudaMatx *J1, *J2;

  /* Linear Stokes I */
  /* Stokes I Mattieu's simplification with h=>X, v=>Y
     | XX XY | =  |gxx*hxx+gxy*hxy  gxy*hyy+gxx*hyx| * v
     | YX YY |    |gyx*hxx+gyy*hxy  gyy*hyy+gyx*hyx|
  cudaMatxZero2C(sumVis);  /* Init output */
  for (i=0; i<numComp; i++) {
    /* Use Mattieu's simplified case with h=>X, v=>Y*/
    /* g??, h?? replaced with pointers into Jones matrices */
    CUDA_COMPLEX_SET(&v,  compFlux[i]*cosArr[i],  compFlux[i]*sinArr[i]) ;/* Visibility model */
    J1 = Jones1[i]; J2 = Jones2[i];
    CUDA_COMPLEX_MULC(&t1, J1->cpx[0], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[1], J2->cpx[1]); /* XX */
    CUDA_COMPLEX_ADD2(&t3, t1, t2);   CUDA_COMPLEX_MUL2(&txx, t3, v);
    CUDA_COMPLEX_MULC(&t1, J1->cpx[1], J2->cpx[3]); CUDA_COMPLEX_MULC(&t2, J1->cpx[0], J2->cpx[2]); /* XY */
    CUDA_COMPLEX_ADD2(&t3, t1, t2);   CUDA_COMPLEX_MUL2(&txy, t3, v);
    CUDA_COMPLEX_MULC(&t1, J1->cpx[2], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[3], J2->cpx[1]); /* YX */
    CUDA_COMPLEX_ADD2(&t3, t1, t2);   CUDA_COMPLEX_MUL2(&tyx, t3, v);
    CUDA_COMPLEX_MULC(&t1, J1->cpx[2], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[3], J2->cpx[3]); /* YY  */
    CUDA_COMPLEX_ADD2(&t3, t1, t2);   CUDA_COMPLEX_MUL2(&tyy, t3, v);
    /* Sum into sumVis */
    CUDA_COMPLEX_ADD2(&sumVis->cpx[0], txx, sumVis->cpx[0]);
    CUDA_COMPLEX_ADD2(&sumVis->cpx[1], txy, sumVis->cpx[1]);
    CUDA_COMPLEX_ADD2(&sumVis->cpx[2], tyx, sumVis->cpx[2]);
    CUDA_COMPLEX_ADD2(&sumVis->cpx[3], tyy, sumVis->cpx[3]);
  } /* end loop over components */
  return; 
} /* End cudaJonesCorSum_Lin_I */

/** Public: Jones correct one SkyModel components, Lin pol feed, Stokes I 
 * Converts SkyModel comps to pseudovisibilities, multiply by Jones matrices and sum.
 * \param flux     Component amplitude with any factors multiplied
 * \param sinPh    Array of sine of phase to apply
 * \param cosPh    Array of cosine of phase to apply
 * \param Jones1   Jones matrix array of first antenna of baseline
 * \param Jones2   Jones matrix array of second antenna of baseline
 * \param c2p      cosine (Parallactic angle*2) Not used for I
 * \param s2p      sine (Parallactic angle*2)       "
 * \param sumVis   [out] 2x2 complex matrix for resultant sum, should be zeroed before 1st call
 */
 extern "C"
__device__ void cudaJonesCor1_Lin_I (float flux, float sinPh, float cosPh,
		                     cudaMatx *Jones1, cudaMatx *Jones2, 
				     float c2p, float s2p,
				     cudaMatx *sumVis)
{
  cudaComplex v, t1, t2, t3, txx, tyy, txy, tyx;
  cudaMatx *J1, *J2;

  /* Linear Stokes I */
  /* Stokes I Mattieu's simplification with h=>X, v=>Y
     | XX XY | =  |gxx*hxx+gxy*hxy  gxy*hyy+gxx*hyx| * v
     | YX YY |    |gyx*hxx+gyy*hxy  gyy*hyy+gyx*hyx|
    /* Use Mattieu's simplified case with h=>X, v=>Y*/
    /* g??, h?? replaced with pointers into Jones matrices */
  CUDA_COMPLEX_SET(&v,  flux*cosPh,  flux*sinPh) ;/* Visibility model */
  J1 = Jones1; J2 = Jones2;
  CUDA_COMPLEX_MULC(&t1, J1->cpx[0], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[1], J2->cpx[1]); /* XX */
  CUDA_COMPLEX_ADD2(&t3, t1, t2);   CUDA_COMPLEX_MUL2(&txx, t3, v);
  CUDA_COMPLEX_MULC(&t1, J1->cpx[1], J2->cpx[3]); CUDA_COMPLEX_MULC(&t2, J1->cpx[0], J2->cpx[2]); /* XY */
  CUDA_COMPLEX_ADD2(&t3, t1, t2);   CUDA_COMPLEX_MUL2(&txy, t3, v);
  CUDA_COMPLEX_MULC(&t1, J1->cpx[2], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[3], J2->cpx[1]); /* YX */
  CUDA_COMPLEX_ADD2(&t3, t1, t2);   CUDA_COMPLEX_MUL2(&tyx, t3, v);
  CUDA_COMPLEX_MULC(&t1, J1->cpx[2], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[3], J2->cpx[3]); /* YY  */
  CUDA_COMPLEX_ADD2(&t3, t1, t2);   CUDA_COMPLEX_MUL2(&tyy, t3, v);
  /* Add into sumVis */
  CUDA_COMPLEX_ADD2(&sumVis->cpx[0], txx, sumVis->cpx[0]);
  CUDA_COMPLEX_ADD2(&sumVis->cpx[1], txy, sumVis->cpx[1]);
  CUDA_COMPLEX_ADD2(&sumVis->cpx[2], tyx, sumVis->cpx[2]);
  CUDA_COMPLEX_ADD2(&sumVis->cpx[3], tyy, sumVis->cpx[3]);
  return; 
} /* End cudaJonesCor1_Lin_I */

/** Public: Jones correct and sum SkyModel components, Lin pol feed, Stokes Q model
 * Converts SkyModel comps to pseudovisibilities, multiply by Jones matrices and sum.
 * \param numComp  Number of SkyModel components
 * \param compFlux Component amplitude with any factors multiplied
 * \param sinArr   Array of sine of phase to apply, per numComp
 * \param cosArr   Array of cosine of phase to apply, per numComp
 * \param Jones1   Jones matrix array of first antenna of baseline, per numComp
 * \param Jones2   Jones matrix array of second antenna of baseline, per numComp
 * \param c2p      cosine (Parallactic angle*2)
 * \param s2p      sine (Parallactic angle*2)
 * \param sumVis   [out] 2x2 complex matrix for resultant sum
 */
 extern "C"
__device__ void cudaJonesCorSum_Lin_Q (long numComp, float *compFlux, 
	                               float *sinArr, float *cosArr,
		                       cudaMatx **Jones1, cudaMatx **Jones2, 
                                       float c2p, float s2p,
			               cudaMatx *sumVis)
{
  long i;
  cudaComplex v, t1, t2, t3, t4, txx, tyy, txy, tyx;
  cudaMatx *J1, *J2;

  /* Linear Stokes Q */
  /* Stokes Q Mattieu's simplification with h=>X, v=>Y
    | XX XY | =
    | YX YY |
      |c2p*(gxx*hxx-gxy*hxy)-s2p*(gxy*hxx+gxx*hxy)  c2p*(gxx*hyx-gxy*hyy)-s2p*(gxy*hyx+gxx*hyy)| * v
      |c2p*(gyx*hxx-gyy*hxy)-s2p*(gyy*hxx+gyx*hxy)  c2p*(gyx*hyx-gyy*hyy)-s2p*(gyy*hyx+gyx*hyy)|      */
  cudaMatxZero2C(sumVis);  /* Init output  */
  for (i=0; i<numComp; i++) {
    /* Use Mattieu's simplified case with h=>X, v=>Y*/
    CUDA_COMPLEX_SET(&v, compFlux[i]*cosArr[i], compFlux[i]*sinArr[i]) ;/* Visibility model */
    J1 = Jones1[i]; J2 = Jones2[i];
    CUDA_COMPLEX_MULC(&t1, J1->cpx[0], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[1], J2->cpx[1]); 
    CUDA_COMPLEX_SUB (&t3, t1, t2); t3.real*=c2p; t3.imag*=c2p; /* XX */
    CUDA_COMPLEX_MULC(&t1, J1->cpx[1], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[0], J2->cpx[1]); 
    CUDA_COMPLEX_ADD2(&t4, t1, t2); t4.real*=s2p; t4.imag*=s2p;
    CUDA_COMPLEX_SUB(&t1, t3, t4); CUDA_COMPLEX_MUL2(&txx, t1, v);
    CUDA_COMPLEX_MULC(&t1, J1->cpx[0], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[1], J2->cpx[3]); 
    CUDA_COMPLEX_SUB (&t3, t1, t2); t3.real*=c2p; t3.imag*=c2p; /* XY */
    CUDA_COMPLEX_MULC(&t1, J1->cpx[1], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[0], J2->cpx[3]); 
    CUDA_COMPLEX_ADD2(&t4, t1, t2); t4.real*=s2p; t4.imag*=s2p;
    CUDA_COMPLEX_SUB(&t1, t3, t4); CUDA_COMPLEX_MUL2(&txy, t1, v);
    CUDA_COMPLEX_MULC(&t1, J1->cpx[2], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[3], J2->cpx[1]); 
    CUDA_COMPLEX_SUB (&t3, t1, t2); t3.real*=c2p; t3.imag*=c2p; /* YX */
    CUDA_COMPLEX_MULC(&t1, J1->cpx[3], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[2], J2->cpx[1]); 
    CUDA_COMPLEX_ADD2(&t4, t1, t2); t4.real*=s2p; t4.imag*=s2p;
    CUDA_COMPLEX_SUB(&t1, t3, t4); CUDA_COMPLEX_MUL2(&tyx, t1, v);
    CUDA_COMPLEX_MULC(&t1, J1->cpx[2], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[3], J2->cpx[3]); 
    CUDA_COMPLEX_SUB (&t3, t1, t2); t3.real*=c2p; t3.imag*=c2p; /* YY */
    CUDA_COMPLEX_MULC(&t1, J1->cpx[3], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[2], J2->cpx[3]); 
    CUDA_COMPLEX_ADD2(&t4, t1, t2); t4.real*=s2p; t4.imag*=s2p;
    CUDA_COMPLEX_SUB(&t1, t3, t4); CUDA_COMPLEX_MUL2(&tyy, t1, v);
    /* Sum into sumVis */
    CUDA_COMPLEX_ADD2(&sumVis->cpx[0], txx, sumVis->cpx[0]);
    CUDA_COMPLEX_ADD2(&sumVis->cpx[1], txy, sumVis->cpx[1]);
    CUDA_COMPLEX_ADD2(&sumVis->cpx[2], tyx, sumVis->cpx[2]);
    CUDA_COMPLEX_ADD2(&sumVis->cpx[3], tyy, sumVis->cpx[3]);
  } /* end loop over components */
  return;  
} /* End cudaJonesCorSum_Lin_Q */

/** Public: Jones correct and sum a single SkyModel component, Lin pol feed, Stokes Q model
 * Converts SkyModel comp to pseudovisibility, multiply by Jones matrices and sum.
 * \param flux     Component amplitude with any factors multiplied
 * \param sinPh    Sine of phase to apply
 * \param cosPh    Cosine of phase to apply
 * \param Jones1   Jones matrix array of first antenna of baseline
 * \param Jones2   Jones matrix array of second antenna of baseline
 * \param c2p      cosine (Parallactic angle*2)
 * \param s2p      sine (Parallactic angle*2)
 * \param sumVis   [out] 2x2 complex matrix for resultant sum, should be zeroed before 1st call
 */
 extern "C"
 __device__ void cudaJonesCor1_Lin_Q (float flux, float sinPh, float cosPh,
		                      cudaMatx *Jones1, cudaMatx *Jones2, 
                                      float c2p, float s2p,
			              cudaMatx *sumVis)
{
  cudaComplex v, t1, t2, t3, t4, txx, tyy, txy, tyx;
  cudaMatx *J1, *J2;

  /* Linear Stokes Q */
  /* Stokes Q Mattieu's simplification with h=>X, v=>Y
    | XX XY | =
    | YX YY |
      |c2p*(gxx*hxx-gxy*hxy)-s2p*(gxy*hxx+gxx*hxy)  c2p*(gxx*hyx-gxy*hyy)-s2p*(gxy*hyx+gxx*hyy)| * v
      |c2p*(gyx*hxx-gyy*hxy)-s2p*(gyy*hxx+gyx*hxy)  c2p*(gyx*hyx-gyy*hyy)-s2p*(gyy*hyx+gyx*hyy)|      */
    /* Use Mattieu's simplified case with h=>X, v=>Y*/
   CUDA_COMPLEX_SET(&v, flux*cosPh, flux*sinPh) ;/* Visibility model */
   J1 = Jones1; J2 = Jones2;
   CUDA_COMPLEX_MULC(&t1, J1->cpx[0], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[1], J2->cpx[1]); CUDA_COMPLEX_SUB (&t3, t1, t2); t3.real*=c2p; t3.imag*=c2p; /* XX */
   CUDA_COMPLEX_MULC(&t1, J1->cpx[1], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[0], J2->cpx[1]); CUDA_COMPLEX_ADD2(&t4, t1, t2); t4.real*=s2p; t4.imag*=s2p;
   CUDA_COMPLEX_SUB(&t1, t3, t4); CUDA_COMPLEX_MUL2(&txx, t1, v);
   
   CUDA_COMPLEX_MULC(&t1, J1->cpx[0], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[1], J2->cpx[3]); CUDA_COMPLEX_SUB (&t3, t1, t2); t3.real*=c2p; t3.imag*=c2p; /* XY */
   CUDA_COMPLEX_MULC(&t1, J1->cpx[1], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[0], J2->cpx[3]); CUDA_COMPLEX_ADD2(&t4, t1, t2); t4.real*=s2p; t4.imag*=s2p;
   CUDA_COMPLEX_SUB(&t1, t3, t4); CUDA_COMPLEX_MUL2(&txy, t1, v);
   
   CUDA_COMPLEX_MULC(&t1, J1->cpx[2], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[3], J2->cpx[1]); CUDA_COMPLEX_SUB (&t3, t1, t2); t3.real*=c2p; t3.imag*=c2p; /* YX */
   CUDA_COMPLEX_MULC(&t1, J1->cpx[3], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[2], J2->cpx[1]); CUDA_COMPLEX_ADD2(&t4, t1, t2); t4.real*=s2p; t4.imag*=s2p;
   CUDA_COMPLEX_SUB(&t1, t3, t4); CUDA_COMPLEX_MUL2(&tyx, t1, v);
   
   CUDA_COMPLEX_MULC(&t1, J1->cpx[2], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[3], J2->cpx[3]); CUDA_COMPLEX_SUB (&t3, t1, t2); t3.real*=c2p; t3.imag*=c2p; /* YY */
   CUDA_COMPLEX_MULC(&t1, J1->cpx[3], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[2], J2->cpx[3]); CUDA_COMPLEX_ADD2(&t4, t1, t2); t4.real*=s2p; t4.imag*=s2p;
   CUDA_COMPLEX_SUB(&t1, t3, t4); CUDA_COMPLEX_MUL2(&tyy, t1, v);
   /* Sum into sumVis */
   CUDA_COMPLEX_ADD2(&sumVis->cpx[0], txx, sumVis->cpx[0]);
   CUDA_COMPLEX_ADD2(&sumVis->cpx[1], txy, sumVis->cpx[1]);
   CUDA_COMPLEX_ADD2(&sumVis->cpx[2], tyx, sumVis->cpx[2]);
   CUDA_COMPLEX_ADD2(&sumVis->cpx[3], tyy, sumVis->cpx[3]);
 return;  
} /* End cudaJonesCor1Sum_Lin_Q */

/** Public: Jones correct and sum SkyModel components, Lin pol feed, Stokes U model
 * converts SkyModel comps to pseudovisibilities, multiply by Jones matrices and sum.
 * \param numComp  Number of SkyModel components
 * \param compFlux Component amplitude with any factors multiplied
 * \param sinArr   Array of sine of phase to apply, per numComp
 * \param cosArr   Array of cosine of phase to apply, per numComp
 * \param Jones1   Jones matrix array of first antenna of baseline, per numComp
 * \param Jones2   Jones matrix array of second antenna of baseline, per numComp
 * \param c2p      cosine (Parallactic angle*2)
 * \param s2p      sine (Parallactic angle*2)
 * \param sumVis   [out] 2x2 complex matrix for resultant sum
 */
 extern "C"
__device__ void cudaJonesCorSum_Lin_U (long numComp, float *compFlux, 
	                              float *sinArr, float *cosArr,
		                      cudaMatx **Jones1, cudaMatx **Jones2, 
				      float c2p, float s2p,
			              cudaMatx *sumVis)
{
  long i;
  cudaComplex v, t1, t2, t3, t4, txx, tyy, txy, tyx;
  cudaMatx *J1, *J2;

  /* Linear Stokes U */
  /* Stokes U Mattieu's simplification with h=>X, v=>Y
     | XX XY | =
     | YX YY |
        |c2p*(gxy*hxx+gxx*hxy)+s2p*(gxx*hxx-gxy*hxy)  c2p*(gxy*hyx+gxx*hyy)+s2p*(gxx*hyx-gxy*hyy)| * v
        |c2p*(gyy*hxx+gyx*hxy)+s2p*(gyx*hxx-gyy*hxy)  c2p*(gyy*hyx+gyx*hyy)+s2p*(gyx*hyx-gyy*hyy)|      */
  cudaMatxZero2C(sumVis);  /* Init output */
  for (i=0; i<numComp; i++) {
    /* Use Mattieu's simplified case with h=>X, v=>Y*/
    CUDA_COMPLEX_SET(&v, compFlux[i]*cosArr[i], compFlux[i]*sinArr[i]) ;/* Visibility */
    J1 = Jones1[i]; J2 = Jones2[i];
    CUDA_COMPLEX_MULC(&t1, J1->cpx[1], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[0], J2->cpx[1]); 
    CUDA_COMPLEX_ADD2(&t3, t1, t2); t3.real*=c2p; t3.imag*=c2p; /* XX */
    CUDA_COMPLEX_MULC(&t1, J1->cpx[0], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[1], J2->cpx[1]); 
    CUDA_COMPLEX_SUB (&t4, t1, t2); t4.real*=s2p; t4.imag*=s2p;
    CUDA_COMPLEX_ADD2(&t1, t3, t4); CUDA_COMPLEX_MUL2(&txx, t1, v);
    CUDA_COMPLEX_MULC(&t1, J1->cpx[1], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[0], J2->cpx[3]); 
    CUDA_COMPLEX_ADD2(&t3, t1, t2); t3.real*=c2p; t3.imag*=c2p; /* XY */
    CUDA_COMPLEX_MULC(&t1, J1->cpx[0], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[1], J2->cpx[3]); 
    CUDA_COMPLEX_SUB (&t4, t1, t2); t4.real*=s2p; t4.imag*=s2p;
    CUDA_COMPLEX_ADD2(&t1, t3, t4); CUDA_COMPLEX_MUL2(&txy, t1, v);
    CUDA_COMPLEX_MULC(&t1, J1->cpx[3], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[2], J2->cpx[1]); 
    CUDA_COMPLEX_ADD2(&t3, t1, t2); t3.real*=c2p; t3.imag*=c2p; /* YX */
    CUDA_COMPLEX_MULC(&t1, J1->cpx[2], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[3], J2->cpx[1]); 
    CUDA_COMPLEX_SUB (&t4, t1, t2); t4.real*=s2p; t4.imag*=s2p;
    CUDA_COMPLEX_ADD2(&t1, t3, t4); CUDA_COMPLEX_MUL2(&tyx, t1, v);
    CUDA_COMPLEX_MULC(&t1, J1->cpx[3], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[2], J2->cpx[3]); 
    CUDA_COMPLEX_ADD2(&t3, t1, t2); t3.real*=c2p; t3.imag*=c2p; /* YY */
    CUDA_COMPLEX_MULC(&t1, J1->cpx[2], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[3], J2->cpx[3]); 
    CUDA_COMPLEX_SUB (&t4, t1, t2); t4.real*=s2p; t4.imag*=s2p;
    CUDA_COMPLEX_ADD2(&t1, t3, t4); CUDA_COMPLEX_MUL2(&tyy, t1, v);
    /* Sum into sumVis */
    CUDA_COMPLEX_ADD2(&sumVis->cpx[0], txx, sumVis->cpx[0]);
    CUDA_COMPLEX_ADD2(&sumVis->cpx[1], txy, sumVis->cpx[1]);
    CUDA_COMPLEX_ADD2(&sumVis->cpx[2], tyx, sumVis->cpx[2]);
    CUDA_COMPLEX_ADD2(&sumVis->cpx[3], tyy, sumVis->cpx[3]);
  } /* end loop over components */
  return;  
} /* End cudaJonesCorSum_Lin_U */

/** Public: Jones correct and sum SkyModel components, Lin pol feed, Stokes U model
 * converts SkyModel comps to pseudovisibilities, multiply by Jones matrices and sum.
 * Converts SkyModel comps to pseudovisibilities, multiply by Jones matrices and sum.
 * \param flux     Component amplitude with any factors multiplied
 * \param sinPh    Array of sine of phase to apply
 * \param cosPh    Array of cosine of phase to apply
 * \param Jones1   Jones matrix array of first antenna of baseline
 * \param Jones2   Jones matrix array of second antenna of baseline
 * \param c2p      cosine (Parallactic angle*2)
 * \param s2p      sine (Parallactic angle*2)
 * \param sumVis   [out] 2x2 complex matrix for resultant sum, should be zeroed before 1st call
 */
 extern "C"
__device__ void cudaJonesCor1_Lin_U (float flux, float sinPh, float cosPh,
		                      cudaMatx *Jones1, cudaMatx *Jones2, 
				      float c2p, float s2p,
			              cudaMatx *sumVis)
{
  cudaComplex v, t1, t2, t3, t4, txx, tyy, txy, tyx;
  cudaMatx *J1, *J2;

  /* Linear Stokes U */
  /* Stokes U Mattieu's simplification with h=>X, v=>Y
     | XX XY | =
     | YX YY |
        |c2p*(gxy*hxx+gxx*hxy)+s2p*(gxx*hxx-gxy*hxy)  c2p*(gxy*hyx+gxx*hyy)+s2p*(gxx*hyx-gxy*hyy)| * v
        |c2p*(gyy*hxx+gyx*hxy)+s2p*(gyx*hxx-gyy*hxy)  c2p*(gyy*hyx+gyx*hyy)+s2p*(gyx*hyx-gyy*hyy)|      */
  /* Use Mattieu's simplified case with h=>X, v=>Y*/
  CUDA_COMPLEX_SET(&v, flux*cosPh, flux*sinPh) ;/* Visibility */
  J1 = Jones1; J2 = Jones2;
  CUDA_COMPLEX_MULC(&t1, J1->cpx[1], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[0], J2->cpx[1]); CUDA_COMPLEX_ADD2(&t3, t1, t2); t3.real*=c2p; t3.imag*=c2p; /* XX */
  CUDA_COMPLEX_MULC(&t1, J1->cpx[0], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[1], J2->cpx[1]); CUDA_COMPLEX_SUB (&t4, t1, t2); t4.real*=s2p; t4.imag*=s2p;
  CUDA_COMPLEX_ADD2(&t1, t3, t4); CUDA_COMPLEX_MUL2(&txx, t1, v);
  CUDA_COMPLEX_MULC(&t1, J1->cpx[1], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[0], J2->cpx[3]); CUDA_COMPLEX_ADD2(&t3, t1, t2); t3.real*=c2p; t3.imag*=c2p; /* XY */
  CUDA_COMPLEX_MULC(&t1, J1->cpx[0], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[1], J2->cpx[3]); CUDA_COMPLEX_SUB (&t4, t1, t2); t4.real*=s2p; t4.imag*=s2p;
  CUDA_COMPLEX_ADD2(&t1, t3, t4); CUDA_COMPLEX_MUL2(&txy, t1, v);
  CUDA_COMPLEX_MULC(&t1, J1->cpx[3], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[2], J2->cpx[1]); CUDA_COMPLEX_ADD2(&t3, t1, t2); t3.real*=c2p; t3.imag*=c2p; /* YX */
  CUDA_COMPLEX_MULC(&t1, J1->cpx[2], J2->cpx[0]); CUDA_COMPLEX_MULC(&t2, J1->cpx[3], J2->cpx[1]); CUDA_COMPLEX_SUB (&t4, t1, t2); t4.real*=s2p; t4.imag*=s2p;
  CUDA_COMPLEX_ADD2(&t1, t3, t4); CUDA_COMPLEX_MUL2(&tyx, t1, v);
  CUDA_COMPLEX_MULC(&t1, J1->cpx[3], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[2], J2->cpx[3]); CUDA_COMPLEX_ADD2(&t3, t1, t2); t3.real*=c2p; t3.imag*=c2p; /* YY */
  CUDA_COMPLEX_MULC(&t1, J1->cpx[2], J2->cpx[2]); CUDA_COMPLEX_MULC(&t2, J1->cpx[3], J2->cpx[3]); CUDA_COMPLEX_SUB (&t4, t1, t2); t4.real*=s2p; t4.imag*=s2p;
  CUDA_COMPLEX_ADD2(&t1, t3, t4); CUDA_COMPLEX_MUL2(&tyy, t1, v);
  /* Sum into sumVis */
  CUDA_COMPLEX_ADD2(&sumVis->cpx[0], txx, sumVis->cpx[0]);
  CUDA_COMPLEX_ADD2(&sumVis->cpx[1], txy, sumVis->cpx[1]);
  CUDA_COMPLEX_ADD2(&sumVis->cpx[2], tyx, sumVis->cpx[2]);
  CUDA_COMPLEX_ADD2(&sumVis->cpx[3], tyy, sumVis->cpx[3]);
  return;  
} /* End cudaJonesCor1_Lin_U */

 /** Public: Jones correct and sum SkyModel components, Lin pol feed, Stokes V model
 * converts SkyModel comps to pseudovisibilities, multiply by Jones matrices and sum.
 * \param numComp  Number of SkyModel components
 * \param compFlux Component amplitude with any factors multiplied
 * \param sinArr   Array of sine of phase to apply, per numComp
 * \param cosArr   Array of cosine of phase to apply, per numComp
 * \param Jones1   Jones matrix array of first antenna of baseline, per numComp
 * \param Jones2   Jones matrix array of second antenna of baseline, per numComp
 * \param c2p      cosine (Parallactic angle*2) not used here
 * \param s2p      sine (Parallactic angle*2) not used here
 * \param sumVis   [out] 2x2 complex matrix for resultant sum
 */
 extern "C"
__device__ void cudaJonesCorSum_Lin_V (long numComp, float *compFlux, 
	                               float *sinArr, float *cosArr,
		                       cudaMatx **Jones1, cudaMatx **Jones2, 
				       float c2p, float s2p,
			               cudaMatx *sumVis)
{
  long i;
  cudaComplex v, t1, t2, t3, J, txx, tyy, txy, tyx;
  cudaMatx *J1, *J2;

  /* Linear Stokes V */
  /* Stokes V Mattieu's simplification with h=>X, v=>Y
      | XX XY | = |(gxx*hxy-gxy*hxx)  (gxx*hyy-gxy*hyx)| * j * v
      | YX YY |   |(gyx*hxy-gyy*hxx)  (gyx*hyy-gyy*hyx)|           */
  cudaMatxZero2C(sumVis);  /* Init output */
  CUDA_COMPLEX_SET(&J, 0.0, 1.0); /* sqrt(-1) */
  for (i=0; i<numComp; i++) {
    /* Use Mattieu's simplified case with h=>X, v=>Y*/
    CUDA_COMPLEX_SET(&v, compFlux[i]*cosArr[i], compFlux[i]*sinArr[i]) ;/* Visibility model */
    J1 = Jones1[i]; J2 = Jones2[i];
    CUDA_COMPLEX_MULC(&t1, J1->cpx[0], J2->cpx[1]); CUDA_COMPLEX_MULC(&t2, J1->cpx[1], J2->cpx[3]); CUDA_COMPLEX_SUB(&t3, t1, t2);  /* XX */
    CUDA_COMPLEX_MUL3(&txx, J, t3, v);
    CUDA_COMPLEX_MULC(&t1, J1->cpx[0], J2->cpx[3]); CUDA_COMPLEX_MULC(&t2, J1->cpx[1], J2->cpx[2]); CUDA_COMPLEX_SUB(&t3, t1, t2);  /* XY */
    CUDA_COMPLEX_MUL3(&txy, J, t3, v);
    CUDA_COMPLEX_MULC(&t1, J1->cpx[2], J2->cpx[1]); CUDA_COMPLEX_MULC(&t2, J1->cpx[3], J2->cpx[0]); CUDA_COMPLEX_SUB(&t3, t1, t2);  /* YX */
    CUDA_COMPLEX_MUL3(&tyx, J, t3, v);
    CUDA_COMPLEX_MULC(&t1, J1->cpx[2], J2->cpx[3]); CUDA_COMPLEX_MULC(&t2, J1->cpx[3], J2->cpx[2]); CUDA_COMPLEX_SUB(&t3, t1, t2);  /* YY */
    CUDA_COMPLEX_MUL3(&tyy, J, t3, v);
    /* Sum into sumVis */
    CUDA_COMPLEX_ADD2(&sumVis->cpx[0], txx, sumVis->cpx[0]);
    CUDA_COMPLEX_ADD2(&sumVis->cpx[1], txy, sumVis->cpx[1]);
    CUDA_COMPLEX_ADD2(&sumVis->cpx[2], tyx, sumVis->cpx[2]);
    CUDA_COMPLEX_ADD2(&sumVis->cpx[3], tyy, sumVis->cpx[3]);
  } /* end loop over components */
    return;  
} /* End cudaJonesCorSum_Lin_V */

/** Public: Jones correct and sum a single SkyModel component, Lin pol feed, Stokes V model
 * converts SkyModel comps to pseudovisibility, multiply by Jones matrices and sum.
 * \param flux     Component amplitude with any factors multiplied
 * \param sinPh    Sine of phase to apply
 * \param cosPh    Cosine of phase to apply
 * \param Jones1   Jones matrix array of first antenna of baseline
 * \param Jones2   Jones matrix array of second antenna of baseline
 * \param c2p      cosine (Parallactic angle*2) not used here
 * \param s2p      sine (Parallactic angle*2) not used here
 * \param sumVis   [out] 2x2 complex matrix for resultant sum, should be zeroed before 1st call
 */
 extern "C"
__device__ void cudaJonesCor1_Lin_V (float flux, float sinPh, float cosPh,
		                     cudaMatx *Jones1, cudaMatx *Jones2, 
				     float c2p, float s2p,
			             cudaMatx *sumVis)
{
  cudaComplex v, t1, t2, t3, J, txx, tyy, txy, tyx;
  cudaMatx *J1, *J2;

  /* Linear Stokes V */
  /* Stokes V Mattieu's simplification with h=>X, v=>Y
      | XX XY | = |(gxx*hxy-gxy*hxx)  (gxx*hyy-gxy*hyx)| * j * v
      | YX YY |   |(gyx*hxy-gyy*hxx)  (gyx*hyy-gyy*hyx)|           */
  CUDA_COMPLEX_SET(&J, 0.0, 1.0); /* sqrt(-1) */
  CUDA_COMPLEX_SET(&v, flux*cosPh, flux*sinPh) ;/* Visibility model */
  J1 = Jones1; J2 = Jones2;
  CUDA_COMPLEX_MULC(&t1, J1->cpx[0], J2->cpx[1]); CUDA_COMPLEX_MULC(&t2, J1->cpx[1], J2->cpx[3]); CUDA_COMPLEX_SUB(&t3, t1, t2);  /* XX */
  CUDA_COMPLEX_MUL3(&txx, J, t3, v);
  CUDA_COMPLEX_MULC(&t1, J1->cpx[0], J2->cpx[3]); CUDA_COMPLEX_MULC(&t2, J1->cpx[1], J2->cpx[2]); CUDA_COMPLEX_SUB(&t3, t1, t2);  /* XY */
  CUDA_COMPLEX_MUL3(&txy, J, t3, v);
  CUDA_COMPLEX_MULC(&t1, J1->cpx[2], J2->cpx[1]); CUDA_COMPLEX_MULC(&t2, J1->cpx[3], J2->cpx[0]); CUDA_COMPLEX_SUB(&t3, t1, t2);  /* YX */
  CUDA_COMPLEX_MUL3(&tyx, J, t3, v);
  CUDA_COMPLEX_MULC(&t1, J1->cpx[2], J2->cpx[3]); CUDA_COMPLEX_MULC(&t2, J1->cpx[3], J2->cpx[2]); CUDA_COMPLEX_SUB(&t3, t1, t2);  /* YY */
  CUDA_COMPLEX_MUL3(&tyy, J, t3, v);
  /* Sum into sumVis */
  CUDA_COMPLEX_ADD2(&sumVis->cpx[0], txx, sumVis->cpx[0]);
  CUDA_COMPLEX_ADD2(&sumVis->cpx[1], txy, sumVis->cpx[1]);
  CUDA_COMPLEX_ADD2(&sumVis->cpx[2], tyx, sumVis->cpx[2]);
  CUDA_COMPLEX_ADD2(&sumVis->cpx[3], tyy, sumVis->cpx[3]);
  return;  
} /* End cudaJonesCor1_Lin_V */

#endif /* CUDALINCORSUM_H */ 
