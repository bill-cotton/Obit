/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2002-2011                                          */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitUV ClassInfo structure     */
/* This is intended to be included in a classInfo structure definition*/
#include "ObitClassDef.h"  /* Parent class ClassInfo definition file */
/** Function pointer to Constructor. */
ObitFArrayCreateFP ObitFArrayCreate;
/** Function pointer to subarray constructor. */
ObitFArraySubArrFP ObitFArraySubArr;
/** Function pointer to transpose constructor. */
ObitFArrayTransposeFP ObitFArrayTranspose;
/** Function pointer to compatability checker. */
ObitFArrayIsCompatableFP ObitFArrayIsCompatable;
/** Function pointer to reallocate structures. */
ObitFArrayReallocFP ObitFArrayRealloc;
/** Function pointer to return index. */
ObitFArrayIndexFP ObitFArrayIndex;
/** Function pointer to find maximum. */
ObitFArrayMaxFP ObitFArrayMax;
/** Function pointer to find maximum abs. value. */
ObitFArrayMaxAbsFP ObitFArrayMaxAbs;
/** Function pointer to find minimum. */
ObitFArrayMinFP ObitFArrayMin;
/** Function pointer to deblank. */
ObitFArrayDeblankFP ObitFArrayDeblank;
/** Function pointer to find  RMS of pixel distribution from histogram.*/
ObitFArrayRMSFP ObitFArrayRMS;
/** Function pointer to find  RMS of pixel distribution.*/
ObitFArrayRawRMSFP ObitFArrayRawRMS;
/** Function pointer to find  RMS about zero*/
ObitFArrayRMS0FP ObitFArrayRMS0;
/** Function pointer to find  RMS with possible quantization */
ObitFArrayRMSQuantFP ObitFArrayRMSQuant;
/** Function pointer to find  Mode of pixel distribution.*/
ObitFArrayModeFP ObitFArrayMode;
/** Function pointer to find Mean of pixel distribution */
ObitFArrayMeanFP ObitFArrayMean;
/** Function pointer to fill elements of an FArray */
ObitFArrayFillFP ObitFArrayFill;
/** Function pointer to negate elements of an FArray */
ObitFArrayNegFP ObitFArrayNeg;
/** Function pointer to sine of elements of an FArray */
ObitFArraySinFP ObitFArraySin;
/** Function pointer to cosine of elements of an FArray */
ObitFArrayCosFP ObitFArrayCos;
/** Function pointer to square root of elements of an FArray */
ObitFArrayCosFP ObitFArraySqrt;
/** Function pointer to sum elements of an FArray */
ObitFArraySumFP ObitFArraySum;
/** Function pointer to number of valid elements in an FArray */
ObitFArrayCountFP ObitFArrayCount;
/** Function pointer to Add a scalar to elements of an FArray */
ObitFArraySAddFP ObitFArraySAdd;
/** Function pointer to Multiply elements of an FArray by a scalar*/
ObitFArraySMulFP ObitFArraySMul;
/** Function pointer to divide elements of an FArray into a scalar*/
ObitFArraySDivFP ObitFArraySDiv;
/** Function pointer to Clip elements of an FArray outside of a given range */
ObitFArrayClipFP ObitFArrayClip;
/** Function pointer to Clip elements of an FArray inside of a given range */
ObitFArrayInClipFP ObitFArrayInClip;
/** Function pointer to Blank elements of an array where another is blanked */
ObitFArrayBlankFP ObitFArrayBlank;
/** Function pointer to get larger elements of two FArrays */
ObitFArrayMaxArrFP ObitFArrayMaxArr;
/** Function pointer to get lesser elements of two FArrays */
ObitFArrayMinArrFP ObitFArrayMinArr;
/** Function pointer to get more extreme elements of two FArrays */
ObitFArrayExtArrFP ObitFArrayExtArr;
/** Function pointer to Add nonblanked elements of two FArrays */
ObitFArraySumArrFP ObitFArraySumArr;
/** Function pointer to average nonblanked elements of two FArrays */
ObitFArrayAvgArrFP ObitFArrayAvgArr;
/** Function pointer to Add elements of two FArrays */
ObitFArrayAddFP ObitFArrayAdd;
/** Function pointer to Subtract elements of two FArrays */
ObitFArraySubFP ObitFArraySub;
/** Function pointer to Multiply elements of two FArrays */
ObitFArrayMulFP ObitFArrayMul;
/** Function pointer to Divide elements of two FArrays */
ObitFArrayDivFP ObitFArrayDiv;
/** Function pointer to Divide elements of two FArrays with clipping */
ObitFArrayDivClipFP ObitFArrayDivClip;
/** Function pointer to "dot" product of two FArrays */
ObitFArrayDotFP ObitFArrayDot;
/** Function pointer to multiple 2D by row, col vector */
ObitFArrayMulColRowFP ObitFArrayMulColRow;
/** Function pointer to determine quantization */
ObitFArrayQuantFP ObitFArrayQuant;
/** Function pointer to  Convert a 1D "center at edges" array to proper order  */
ObitFArray1DCenterFP ObitFArray1DCenter;
/** Function pointer to  Convert a 2D "center at edges" array to proper order  */
ObitFArray2DCenterFP ObitFArray2DCenter;
/** Function pointer to  inplace invert a symmetric 2D array  */
ObitFArray2DSymInvFP ObitFArray2DSymInv;
/** Function pointer to  Make 2-D Circular Gaussian  */
ObitFArray2DCGaussFP ObitFArray2DCGauss;
/** Function pointer to  Make 2-D Elliptical Gaussian  */
ObitFArray2DCGaussFP ObitFArray2DEGauss;
/** Function pointer to Shift and Add scaled array  */
ObitFArrayShiftAddFP ObitFArrayShiftAdd;
/** Function pointer to Zero pad an array  */
ObitFArrayPadFP ObitFArrayPad;
/** Function pointer to Convolve a list of Gaussians onto an FArray   */
ObitFArrayConvGausFP ObitFArrayConvGaus;
/** Function pointer to Select elements in an FArray by increment */
ObitFArraySelIncFP ObitFArraySelInc;
/** Function pointer to return histogram of elements in an FArray */
ObitFArrayHistoFP ObitFArrayHisto;
