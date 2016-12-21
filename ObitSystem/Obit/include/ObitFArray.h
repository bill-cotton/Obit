/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2016                                          */
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
#ifndef OBITFARRAY_H 
#define OBITFARRAY_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitInfoList.h"
#include "ObitThread.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFArray.h
 * ObitFArray numeric array class definition.
 *
 * This class is derived from the #Obit class.
 * Related functions are in the 
 * \link ObitFArrayUtil.h ObitFArrayUtil 
 * \endlink module.
 *
 * This class is for creating and manipulating a Array as a memory resident 
 * multidimensional rectangular array of floats.
 * Elements are stored in order of the increasing axis order (the reverse of the
 * usual c definition).
 * Except as noted, magic value blanking is supported (OBIT_MAGIC).
 * 
 * \section ObitFArrayaccess Creators and Destructors
 * An ObitFArray will usually be created using ObitFArrayCreate which allows 
 * specifying a name for the object as well as dimensionality of the array.
 *
 * A copy of a pointer to an ObitFArray should always be made using the
 * #ObitFArrayRef function which updates the reference count in the object.
 * Then whenever freeing an ObitFArray or changing a pointer, the function
 * #ObitFArrayUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitFArray Class structure. */
typedef struct {
#include "ObitFArrayDef.h"   /* this class definition */
} ObitFArray;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitFArray
 * returns a ObitFArray*.
 * in = object to unreference
 */
#define ObitFArrayUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitFArray.
 * returns a ObitFArray*.
 * in = object to reference
 */
#define ObitFArrayRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitFArrayIsA(in) ObitIsA (in, ObitFArrayGetClass())

/** Maximum ObitFArray number of dimensions */
#ifndef MAXFARRAYDIM
#define MAXFARRAYDIM 10
#endif

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitFArrayClassInit (void);

/** Public: Default Constructor. */
ObitFArray* newObitFArray (gchar* name);

/** Public: Create/initialize ObitFArray structures */
ObitFArray* ObitFArrayCreate (gchar* name, olong ndim, olong *naxis);
/** Typedef for definition of class pointer structure */
typedef void (*ObitFArrayCreateFP) (gchar* name, olong ndim, olong *naxis);

/** Public: ClassInfo pointer */
gconstpointer ObitFArrayGetClass (void);

/** Public: Copy (deep) constructor. */
ObitFArray* ObitFArrayCopy  (ObitFArray *in, ObitFArray *out, ObitErr *err);

/** Public: Copy structure. */
void ObitFArrayClone (ObitFArray *in, ObitFArray *out, ObitErr *err);

/** Public: Copy Subarray constructor. */
ObitFArray* ObitFArraySubArr  (ObitFArray *in, olong *blc, olong *trc, 
			       ObitErr *err);
typedef ObitFArray* (*ObitFArraySubArrFP) (ObitFArray *in, olong *blc, olong *trc, 
					   ObitErr *err);

/** Public: Transpose constructor. */
ObitFArray* ObitFArrayTranspose  (ObitFArray *in, olong *order, ObitErr *err);
typedef ObitFArray* (*ObitFArrayTransposeFP) (ObitFArray *in, olong *order, 
					      ObitErr *err);

/** Public: Are two FArrays of compatable geometry. */
gboolean ObitFArrayIsCompatable  (ObitFArray *in1, ObitFArray *in2);
typedef gboolean (*ObitFArrayIsCompatableFP) (ObitFArray *in1, ObitFArray *in2);

/** Public: Reallocate/initialize ObitFArray structures */
ObitFArray* ObitFArrayRealloc (ObitFArray* in, olong ndim, olong *naxis);
typedef void (*ObitFArrayReallocFP) (ObitFArray* in, olong ndim, olong *naxis);

/** Public: return pointer to a specified element */
ofloat* ObitFArrayIndex (ObitFArray* in, olong *pos);
typedef ofloat* (*ObitFArrayIndexFP) (ObitFArray* in, olong *pos);

/** Public: Find Maximum value in an ObitFArray */
ofloat ObitFArrayMax (ObitFArray* in, olong *pos);
typedef ofloat (*ObitFArrayMaxFP) (ObitFArray* in, olong *pos);

/** Public: Find Maximum abs value in an ObitFArray */
ofloat ObitFArrayMaxAbs (ObitFArray* in, olong *pos);
typedef ofloat (*ObitFArrayMaxAbsFP) (ObitFArray* in, olong *pos);

/** Public: Find Minimum value in an ObitFArray */
ofloat ObitFArrayMin (ObitFArray* in, olong *pos);
typedef ofloat (*ObitFArrayMinFP) (ObitFArray* in, olong *pos);

/** Public: Replace blanks in an ObitFArray */
void ObitFArrayDeblank (ObitFArray* in, ofloat scalar);
typedef void (*ObitFArrayDeblankFP) (ObitFArray* in, ofloat scalar);

/** Public: RMS of pixel distribution from histogram */
ofloat ObitFArrayRMS (ObitFArray* in);
typedef ofloat (*ObitFArrayRMSFP) (ObitFArray* in);

/** Public: RMS of pixel distribution. */
ofloat ObitFArrayRawRMS (ObitFArray* in);
typedef ofloat (*ObitFArrayRawRMSFP) (ObitFArray* in);

/** Public: RMS of pixel about zero. */
ofloat ObitFArrayRMS0 (ObitFArray* in);
typedef ofloat (*ObitFArrayRMS0FP) (ObitFArray* in);

/** Public: Function for combining planes of multi frequency images */
ofloat ObitFArrayRMS0ST (ObitFArray* in);
typedef ofloat (*ObitFArrayRMS0STFP) (ObitFArray* in);

/** Public: RMS of pixel in potentially quantized image. */
ofloat ObitFArrayRMSQuant (ObitFArray* in);
typedef ofloat (*ObitFArrayRMSQuantFP) (ObitFArray* in);

/** Public: Determine quantization and offset in an image */
void ObitFArrayQuant (ObitFArray* in, ofloat *quant, ofloat *zero);
typedef void (*ObitFArrayQuantFP) (ObitFArray* in, ofloat *quant, ofloat *zero);

/** Public: Mode of pixel distribution. */
ofloat ObitFArrayMode (ObitFArray* in);
typedef ofloat (*ObitFArrayModeFP) (ObitFArray* in);

/** Public: Mean of pixel distribution. */
ofloat ObitFArrayMean (ObitFArray* in);
typedef ofloat (*ObitFArrayMeanFP) (ObitFArray* in);

/** Public: fill elements of an FArray */
void ObitFArrayFill (ObitFArray* in, ofloat scalar);
typedef void (*ObitFArrayFillFP) (ObitFArray* in, ofloat scalar);

/** Public: negate elements of an FArray */
void ObitFArrayNeg (ObitFArray* in);
typedef void (*ObitFArrayNegFP) (ObitFArray* in);

/** Public: taks absolute value of  elements of an FArray */
void ObitFArrayAbs (ObitFArray* in);
typedef void (*ObitFArrayAbsFP) (ObitFArray* in);

/** Public: sine of  elements of an FArray */
void ObitFArraySin (ObitFArray* in);
typedef void (*ObitFArraySinFP) (ObitFArray* in);

/** Public: cosine of  elements of an FArray */
void ObitFArrayCos (ObitFArray* in);
typedef void (*ObitFArrayCosFP) (ObitFArray* in);

/** Public: square root of  elements of an FArray */
void ObitFArraySqrt (ObitFArray* in);
typedef void (*ObitFArraySqrtFP) (ObitFArray* in);

/** Public: sum elements of an FArray */
ofloat ObitFArraySum (ObitFArray* in);
typedef ofloat (*ObitFArraySumFP) (ObitFArray* in);

/** Public: number of valid elements in an FArray */
olong ObitFArrayCount (ObitFArray* in);
typedef olong (*ObitFArrayCountFP) (ObitFArray* in);

/** Public: Add a scalar to elements of an FArray */
void ObitFArraySAdd (ObitFArray* in, ofloat scalar);
typedef void (*ObitFArraySAddFP) (ObitFArray* in, ofloat scalar);

/** Public: Multiply elements of an FArray by a scalar*/
void ObitFArraySMul (ObitFArray* in, ofloat scalar);
typedef void (*ObitFArraySMulFP) (ObitFArray* in, ofloat scalar);

/** Public: Divide elements of an FArray into a scalar*/
void ObitFArraySDiv (ObitFArray* in, ofloat scalar);
typedef void (*ObitFArraySDivFP) (ObitFArray* in, ofloat scalar);

/** Public: Clip elements of an FArray outside of a given range */
void ObitFArrayClip (ObitFArray* in, ofloat minVal, ofloat maxVal, ofloat newVal);
typedef void (*ObitFArrayClipFP) (ObitFArray* in, ofloat minVal, ofloat maxVal, 
				  ofloat newVal);

/** Public: Clip elements of an FArray inside of a given range */
void ObitFArrayInClip (ObitFArray* in, ofloat minVal, ofloat maxVal, ofloat newVal);
typedef void (*ObitFArrayInClipFP) (ObitFArray* in, ofloat minVal, ofloat maxVal, 
				    ofloat newVal);

/** Public: Blank elements of an array where another is blanked */
void ObitFArrayBlank (ObitFArray* in1, ObitFArray* in2, ObitFArray* out);
typedef void (*ObitFArrayBlankFP) (ObitFArray* in1, ObitFArray* in2, 
				  ObitFArray* out);

/** Public: Get larger elements of two FArrays */
void ObitFArrayMaxArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out);
typedef void (*ObitFArrayMaxArrFP) (ObitFArray* in1, ObitFArray* in2, 
				    ObitFArray* out);

/** Public: Get lesser elements of two FArrays */
void ObitFArrayMinArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out);
typedef void (*ObitFArrayMinArrFP) (ObitFArray* in1, ObitFArray* in2, 
				    ObitFArray* out);

/** Public: Get more extreme elements of two FArrays */
void ObitFArrayExtArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out);
typedef void (*ObitFArrayExtArrFP) (ObitFArray* in1, ObitFArray* in2, 
				    ObitFArray* out);

/** Public: Sum nonblanked elements of two FArrays */
void ObitFArraySumArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out);
typedef void (*ObitFArraySumArrFP) (ObitFArray* in1, ObitFArray* in2, 
				    ObitFArray* out);

/** Public: Average nonblanked elements of two FArrays */
void ObitFArrayAvgArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out);
typedef void (*ObitFArrayAvgArrFP) (ObitFArray* in1, ObitFArray* in2, 
				    ObitFArray* out);

/** Public: Add elements of two FArrays */
void ObitFArrayAdd (ObitFArray* in1, ObitFArray* in2, ObitFArray* out);
typedef void (*ObitFArrayAddFP) (ObitFArray* in1, ObitFArray* in2, 
				  ObitFArray* out);

/** Public: Abs Add elements of two FArrays */
void ObitFArrayAddAbs (ObitFArray* in1, ObitFArray* in2, ObitFArray* out);
typedef void (*ObitFArrayAddAbsFP) (ObitFArray* in1, ObitFArray* in2, 
			  	  ObitFArray* out);

/** Public: Subtract elements of two FArrays */
void ObitFArraySub (ObitFArray* in1, ObitFArray* in2, ObitFArray* out);
typedef void (*ObitFArraySubFP) (ObitFArray* in1, ObitFArray* in2, 
				  ObitFArray* out);

/** Public: Give the elements of one array the sign of the other */
void ObitFArraySign (ObitFArray* in1, ObitFArray* in2);
typedef void (*ObitFArraySignFP) (ObitFArray* in1, ObitFArray* in2);

/** Public: Multiply elements of two FArrays */
void ObitFArrayMul (ObitFArray* in1, ObitFArray* in2, ObitFArray* out);
typedef void (*ObitFArrayMulFP) (ObitFArray* in1, ObitFArray* in2, 
				  ObitFArray* out);

/** Public: Divide elements of two FArrays */
void ObitFArrayDiv (ObitFArray* in1, ObitFArray* in2, ObitFArray* out);
typedef void (*ObitFArrayDivFP) (ObitFArray* in1, ObitFArray* in2, 
				  ObitFArray* out);

/** Public: Divide elements of two FArrays with clipping*/
void ObitFArrayDivClip (ObitFArray* in1, ObitFArray* in2, ofloat minVal, ObitFArray* out);
typedef void (*ObitFArrayDivClipFP) (ObitFArray* in1, ObitFArray* in2, 
				     ofloat minVal, ObitFArray* out);

/** Public: "Dot" product to two arrays */
ofloat ObitFArrayDot (ObitFArray* in1, ObitFArray* in2);
typedef ofloat (*ObitFArrayDotFP) (ObitFArray* in1, ObitFArray* in2);

/** Public: Multiply a 2D array by a Col vector * Row vector */
void ObitFArrayMulColRow (ObitFArray* in, ObitFArray* row, ObitFArray* col,
			  ObitFArray* out);
typedef void (*ObitFArrayMulColRowFP) (ObitFArray* in, ObitFArray* row, 
				       ObitFArray* col, ObitFArray* out);

/** Public: Convert a 1D "center at edges" array to proper order */
void ObitFArray1DCenter (ObitFArray* in);
typedef void (*ObitFArray1DCenterFP) (ObitFArray* in);

/** Public: Convert a 2D "center at edges" array to proper order */
void ObitFArray2DCenter (ObitFArray* in);
typedef void (*ObitFArray2DCenterFP) (ObitFArray* in);

/** Public: inplace invert a symmetric 2D array */
void ObitFArray2DSymInv (ObitFArray* in, olong *ierr);
typedef void (*ObitFArray2DSymInvFP) (ObitFArray* in, olong *ierr);

/** Public: Make 2-D Circular Gaussian in FArray */
void ObitFArray2DCGauss (ObitFArray* in, olong Cen[2], ofloat FWHM);
typedef void (*ObitFArray2DCGaussFP) (ObitFArray* in, olong Cen[2], ofloat FWHM);

/** Public: Make 2-D Eliptical Gaussian in FArray */
void ObitFArray2DEGauss (ObitFArray* in, ofloat amp, ofloat Cen[2], ofloat GauMod[3]);
typedef void (*ObitFArray2DEGaussFP) (ObitFArray* in, ofloat amp, ofloat Cen[2], 
				      ofloat GauMod[3] );

/** Public: Shift and Add scaled array */
void ObitFArrayShiftAdd (ObitFArray* in1, olong *pos1, 
			 ObitFArray* in2, olong *pos2, 
			 ofloat scalar, ObitFArray* out);
typedef void 
(*ObitFArrayShiftAddFP) (ObitFArray* in1, olong *pos1, 
			 ObitFArray* in2, olong *pos2, 
			 ofloat scalar, ObitFArray* out);

/** Public: Zero pad an array */
void  ObitFArrayPad (ObitFArray* in, ObitFArray* out, ofloat factor);
typedef void (*ObitFArrayPadFP) (ObitFArray* in, ObitFArray* out, 
				 ofloat factor);

/** Public: Convolve a list of Gaussians onto an FArray */
void  ObitFArrayConvGaus (ObitFArray* in, ObitFArray* list, olong ncomp, 
			  ofloat gauss[3]);
typedef void (*ObitFArrayConvGausFP) (ObitFArray* in, ObitFArray* list, 
				      olong ncomp, ofloat gauss[3]);

/** Public: Select elements in an FArray by increment */
void  ObitFArraySelInc (ObitFArray* in, ObitFArray* out, olong *blc, olong *trc, 
			olong* inc, ObitErr *err);
typedef void (*ObitFArraySelIncFP) (ObitFArray* in, ObitFArray* out, 
				    olong *blc, olong *trc, olong* inc, ObitErr *err);

/** Public: return histogram of elements in an FArray */
ObitFArray*  ObitFArrayHisto (ObitFArray* in, olong n, ofloat min, ofloat max);
typedef ObitFArray*  (*ObitFArrayHistoFP) (ObitFArray* in, olong n, ofloat min, ofloat max);

/** Public: exponentiate elements in an FArray */
void ObitFArrayExp (ObitFArray* in, ObitFArray* out);
typedef void  (*ObitFArrayExpFP) (ObitFArray* in, ObitFArray* out);

/** Public: natural log of elements in an FArray */
void ObitFArrayLog (ObitFArray* in, ObitFArray* out);
typedef void  (*ObitFArrayLogFP) (ObitFArray* in, ObitFArray* out);

/** Public: natural log of elements in an FArray */
void ObitFArrayPow (ObitFArray* in1, ObitFArray* in2, ObitFArray* out);
typedef void  (*ObitFArrayPowFP) (ObitFArray* in1, ObitFArray* in2, ObitFArray* out);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitFArrayClassDef.h"
} ObitFArrayClassInfo; 

#endif /* OBITFARRAY_H */ 
