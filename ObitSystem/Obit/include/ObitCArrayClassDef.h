/* $Id: ObitCArrayClassDef.h,v 1.4 2005/11/09 12:53:29 bcotton Exp $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2002-2008                                          */
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
/*  Define the basic components of the ObitUV ClassInfo structure     */
/* This is intended to be included in a classInfo structure definition*/
#include "ObitClassDef.h"  /* Parent class ClassInfo definition file */
/** Function pointer to Constructor. */
ObitCArrayCreateFP ObitCArrayCreate;
/** Function pointer to compatability checker. */
ObitCArrayIsCompatableFP ObitCArrayIsCompatable;
/** Function pointer to reallocate structures. */
ObitCArrayReallocFP ObitCArrayRealloc;
/** Function pointer to return index. */
ObitCArrayIndexFP ObitCArrayIndex;
/** Function pointer to find max abs value of a CArray */
ObitCArrayMaxAbsFP ObitCArrayMaxAbs;
/** Function pointer to find Minimum real or imaginary of a CArray */
ObitCArrayMinFP ObitCArrayMin;
/** Function pointer to negate elements of a CArray*/
ObitCArrayNegFP ObitCArrayNeg;
/** Function pointer to complex conjugate elements of a CArray */
ObitCArrayConjgFP ObitCArrayConjg;
/** Function pointer to Fill a CArray with a complex scalar */
ObitCArrayFillFP ObitCArrayFill;
/** Function pointer to Add a scalar to elements of a CArray */
ObitCArraySAddFP ObitCArraySAdd;
/** Function pointer to Multiply elements of a CArray by a scalar*/
ObitCArraySMulFP ObitCArraySMul;
/** Function pointer to Add a complex scalar to elements of a CArray */
ObitCArrayCSAddFP ObitCArrayCSAdd;
/** Function pointer to Multiply elements of a CArray by a complex sscalar*/
ObitCArrayCSMulFP ObitCArrayCSMul;
/** Function pointer to Add elements of two CArrays */
ObitCArrayAddFP ObitCArrayAdd;
/** Function pointer to Subtract elements of two CArrays */
ObitCArraySubFP ObitCArraySub;
/** Function pointer to Multiply elements of two CArrays */
ObitCArrayMulFP ObitCArrayMul;
/** Function pointer to Divide elements of two CArrays */
ObitCArrayDivFP ObitCArrayDiv;

/* CArray - FArray functions */
/** Public: Create an FArray with the same geometry as a CArray. */
ObitCArrayMakeFFP ObitCArrayMakeF;
/** Public: Create an CArray with the same geometry as a FArray. */
ObitCArrayMakeCFP ObitCArrayMakeC;
/** Public: Are an FArray and a CArray of compatable geometry. */
ObitCArrayIsFCompatableFP ObitCArrayIsFCompatable;
/** Public: Multiply a CArray by an FArray */
ObitCArrayFMulFP ObitCArrayFMul;
/** Public: Divide a CArray by an FArray */
ObitCArrayFDivFP ObitCArrayFDiv;
/** Public: Rotate a CArray by an FArray */
ObitCArrayFRotFP ObitCArrayFRot;
/** Public: Add an FArray to a CArray (real part) */
ObitCArrayFAddFP ObitCArrayFAdd;
/** Public: Form A CArray from two FArrays */
ObitCArrayComplexFP ObitCArrayComplex;
/** Public: Return the real elements of a CArray in an FArray */
ObitCArrayRealFP ObitCArrayReal;
/** Public: Return the imaginary elements of a CArray in an FArray */
ObitCArrayImagFP ObitCArrayImag;
/** Public: Return the amplitudes of a CArray in an FArray */
ObitCArrayAmpFP ObitCArrayAmp;
/** Public: Return the phases of a CArray in an FArray */
ObitCArrayPhaseFP ObitCArrayPhase;
/** Function pointer to  Convert a 2D "center at edges" array to proper order  */
ObitCArray2DCenterFP ObitCArray2DCenter;
/** Function pointer to Add conjugate rows to half plane complex image  */
ObitCArrayAddConjgFP ObitCArrayAddConjg;
