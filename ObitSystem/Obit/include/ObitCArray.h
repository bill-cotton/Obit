/* $Id$          */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITCARRAY_H 
#define OBITCARRAY_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitFArray.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCArray.h
 * ObitCArray complex numeric array class definition.
 *
 * This class is derived from the #Obit class.
 *
 * This class is for creating and manipulating a Array as a memory resident 
 * multidimensional rectangular array of complex (pair of floats as (real,imag).
 * Elements are stored in order of the increasing axis order (the reverse of the
 * usual c definition).
 * 
 * \section ObitCArrayaccess Creators and Destructors
 * An ObitCArray will usually be created using ObitCArrayCreate which allows 
 * specifying a name for the object as well as dimensionality of the array.
 *
 * A copy of a pointer to an ObitCArray should always be made using the
 * #ObitCArrayRef function which updates the reference count in the object.
 * Then whenever freeing an ObitCArray or changing a pointer, the function
 * #ObitCArrayUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitCArray Class structure. */
typedef struct {
#include "ObitCArrayDef.h"   /* this class definition */
} ObitCArray;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitCArray
 * returns a ObitCArray*.
 * in = object to unreference
 */
#define ObitCArrayUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitCArray.
 * returns a ObitCArray*.
 * in = object to reference
 */
#define ObitCArrayRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitCArrayIsA(in) ObitIsA (in, ObitCArrayGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitCArrayClassInit (void);

/** Public: Default Constructor. */
ObitCArray* newObitCArray (gchar* name);

/** Public: Create/initialize ObitCArray structures */
ObitCArray* ObitCArrayCreate (gchar* name, olong ndim, olong *naxis);
/** Typedef for definition of class pointer structure */
typedef void (*ObitCArrayCreateFP) (gchar* name, olong ndim, olong *naxis);

/** Public: ClassInfo pointer */
gconstpointer ObitCArrayGetClass (void);

/** Public: Copy (deep) constructor. */
ObitCArray* ObitCArrayCopy  (ObitCArray *in, ObitCArray *out, ObitErr *err);

/** Public: Are two CArrays of compatable geometry. */
gboolean ObitCArrayIsCompatable  (ObitCArray *in1, ObitCArray *in2);
typedef gboolean (*ObitCArrayIsCompatableFP) (ObitCArray *in1, ObitCArray *in2);

/** Public: Reallocate/initialize ObitCArray structures */
ObitCArray* ObitCArrayRealloc (ObitCArray* in, olong ndim, olong *naxis);
typedef void (*ObitCArrayReallocFP) (ObitCArray* in, olong ndim, olong *naxis);

/** Public: return pointer to a specified element */
ofloat* ObitCArrayIndex (ObitCArray* in, olong *pos);
typedef ofloat* (*ObitCArrayIndexFP) (ObitCArray* in, olong *pos);

/** Public: Find Maximum abs value in an ObitCArray */
ofloat ObitCArrayMaxAbs (ObitCArray* in, olong *pos);
typedef ofloat (*ObitCArrayMaxAbsFP) (ObitCArray* in, olong *pos);

/** Public: Find Minimum real or imaginary in an ObitCArray */
ofloat ObitCArrayMin (ObitCArray* in, olong *pos);
typedef ofloat (*ObitCArrayMinFP) (ObitCArray* in, olong *pos);

/** Public: negate elements of an CArray */
void ObitCArrayNeg (ObitCArray* in);
typedef void (*ObitCArrayNegFP) (ObitCArray* in);

/** Public: conjugate elements of an CArray */
void ObitCArrayConjg (ObitCArray* in);
typedef void (*ObitCArrayConjgFP) (ObitCArray* in);

/** Public: Fill a CArray with a complex scalar */
void ObitCArrayFill (ObitCArray* in, ofloat cmpx[2]);
typedef void (*ObitCArrayFillFP) (ObitCArray* in, ofloat cmpx[2]);

/** Public: Add a scalar to elements of a CArray */
void ObitCArraySAdd (ObitCArray* in, ofloat scalar);
typedef void (*ObitCArraySAddFP) (ObitCArray* in, ofloat scalar);

/** Public: Multiply elements of a CArray by a scalar*/
void ObitCArraySMul (ObitCArray* in, ofloat scalar);
typedef void (*ObitCArraySMulFP) (ObitCArray* in, ofloat scalar);

/** Public: Add a complex scalar to elements of a CArray */
void ObitCArrayCSAdd (ObitCArray* in, ofloat scalar[2]);
typedef void (*ObitCArrayCSAddFP) (ObitCArray* in, ofloat scalar[2]);

/** Public: Multiply elements of a CArray by a complex scalar*/
void ObitCArrayCSMul (ObitCArray* in, ofloat scalar[2]);
typedef void (*ObitCArrayCSMulFP) (ObitCArray* in, ofloat scalar[2]);

/** Public: Add elements of two CArrays */
void ObitCArrayAdd (ObitCArray* in1, ObitCArray* in2, ObitCArray* out);
typedef void (*ObitCArrayAddFP) (ObitCArray* in1, ObitCArray* in2, 
				  ObitCArray* out);

/** Public: Subtract elements of two CArrays */
void ObitCArraySub (ObitCArray* in1, ObitCArray* in2, ObitCArray* out);
typedef void (*ObitCArraySubFP) (ObitCArray* in1, ObitCArray* in2, 
				  ObitCArray* out);

/** Public: Multiply elements of two CArrays */
void ObitCArrayMul (ObitCArray* in1, ObitCArray* in2, ObitCArray* out);
typedef void (*ObitCArrayMulFP) (ObitCArray* in1, ObitCArray* in2, 
				  ObitCArray* out);

/** Public: Divide elements of two CArrays */
void ObitCArrayDiv (ObitCArray* in1, ObitCArray* in2, ObitCArray* out);
typedef void (*ObitCArrayDivFP) (ObitCArray* in1, ObitCArray* in2, 
				  ObitCArray* out);

/* CArray - FArray functions */

/** Public: Create an FArray with the same geometry as a CArray. */
ObitFArray* ObitCArrayMakeF  (ObitCArray *Cin);
typedef ObitFArray* (*ObitCArrayMakeFFP) (ObitCArray *Cin);

/** Public: Create a CArray with the same geometry as a FArray. */
ObitCArray* ObitCArrayMakeC  (ObitFArray *Fin);
typedef ObitCArray* (*ObitCArrayMakeCFP) (ObitFArray *Fin);

/** Public: Are an FArray and a CArray of compatable geometry. */
gboolean ObitCArrayIsFCompatable  (ObitCArray *Cin, ObitFArray *Fin);
typedef gboolean (*ObitCArrayIsFCompatableFP) (ObitCArray *Cin, ObitFArray *Fin);

/** Public: Multiply a CArray by an FArray */
void ObitCArrayFMul (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out);
typedef void (*ObitCArrayFMulFP) (ObitCArray* Cin, ObitFArray* Fin, 
				  ObitCArray* out);

/** Public: Divide a CArray by an FArray */
void ObitCArrayFDiv (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out);
typedef void (*ObitCArrayFDivFP) (ObitCArray* Cin, ObitFArray* Fin, 
				  ObitCArray* out);

/** Public: Rotate phases of a CArray ob phases in an FArray */
void ObitCArrayFRot (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out);
typedef void (*ObitCArrayFRotFP) (ObitCArray* Cin, ObitFArray* Fin, 
				  ObitCArray* out);

/** Public: Add an FArray to a CArray (real part) */
void ObitCArrayFAdd (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out);
typedef void (*ObitCArrayFAddFP) (ObitCArray* Cin, ObitFArray* Fin, 
				  ObitCArray* out);

/** Public: Form A CArray from two FArrays */
void ObitCArrayComplex (ObitFArray* Fin1, ObitFArray* Fin2, ObitCArray* out);
typedef void (*ObitCArrayComplexFP) (ObitFArray* Fin1, ObitFArray* Fin2, 
				      ObitCArray* out);

/** Public: Return the real elements of a CArray in an FArray */
void ObitCArrayReal (ObitCArray* in, ObitFArray* out);
typedef void (*ObitCArrayRealFP) (ObitCArray* in, ObitCArray* out);

/** Public: Return the imaginary elements of a CArray in an FArray */
void ObitCArrayImag (ObitCArray* in, ObitFArray* out);
typedef void (*ObitCArrayImagFP) (ObitCArray* in, ObitCArray* out);

/** Public: Return the amplitudes of a CArray in an FArray */
void ObitCArrayAmp (ObitCArray* in, ObitFArray* out);
typedef void (*ObitCArrayAmpFP) (ObitCArray* in, ObitCArray* out);

/** Public: Return the phasess of a CArray in an FArray */
void ObitCArrayPhase (ObitCArray* in, ObitFArray* out);
typedef void (*ObitCArrayPhaseFP) (ObitCArray* in, ObitCArray* out);

/** Public: Convert a 2D "center at edges" array to proper order */
void ObitCArray2DCenter (ObitCArray* in);
typedef void (*ObitCArray2DCenterFP) (ObitCArray* in);


/** Public: Add conjugate columns to half plane complex image */
ObitCArray* ObitCArrayAddConjg (ObitCArray* in, olong numConjCol);
typedef ObitCArray* (*ObitCArrayAddConjgFP) (ObitCArray* in, 
					     olong numConjCol);


/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitCArrayClassDef.h"
} ObitCArrayClassInfo; 

#endif /* OBITCARRAY_H */ 
