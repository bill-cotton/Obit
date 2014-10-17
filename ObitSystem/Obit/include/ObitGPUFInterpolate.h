/* $Id:  $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014                                               */
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
#ifndef OBITGPUSKYMODEL_H 
#define OBITGPUSKYMODEL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitUV.h"
#include "CUDAFInterpolate.h"
#include "ObitFArray.h"
#include "ObitGPUFArray.h"
#define STREAM_COUNT 4  /* Number of GPU Streams */
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGPUFInterpolate.h
 *
 * ObitGPUFInterpolate GPU enhanced FArray interpolation class
 * Uses functions in ObitCUDAFInterpolate
 *
 * 2D FArrays are interpolated using a GPU.
 * 
 * \section ObitGPUFInterpolateaccess Creators and Destructors
 * An ObitGPUFInterpolate will usually be created using ObitGPUFInterpolateCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitGPUFInterpolate should always be made using the
 * #ObitGPUFInterpolateRef function which updates the reference count in the object.
 * Then whenever freeing an ObitGPUFInterpolate or changing a pointer, the function
 * #ObitGPUFInterpolateUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitGPUFInterpolate Class structure. */
typedef struct {
#include "ObitGPUFInterpolateDef.h"   /* this class definition */
} ObitGPUFInterpolate;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitGPUFInterpolate
 * returns a ObitGPUFInterpolate*.
 * in = object to unreference
 */
#define ObitGPUFInterpolateUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitGPUFInterpolate.
 * returns a ObitGPUFInterpolate*.
 * in = object to reference
 */
#define ObitGPUFInterpolateRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGPUFInterpolateIsA(in) ObitIsA (in, ObitGPUFInterpolateGetClass())

/*---------------------------------- Structures ----------------------------*/
/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitGPUFInterpolateClassInit (void);

/** Public: Default Constructor. */
ObitGPUFInterpolate* newObitGPUFInterpolate (gchar* name);

/** Public: Create/initialize ObitGPUFInterpolate structures */
ObitGPUFInterpolate* ObitGPUFInterpolateCreate (gchar* name, ObitFArray *inArray, 
						ObitFArray *xArray, ObitFArray *yArray,
						int hwidth, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitGPUFInterpolate* 
(*ObitGPUFInterpolateCreateFP) (gchar* name, ObitFArray *inArray, 
				ObitFArray *xArray, ObitFArray *yArray,
				int hwidth, ObitErr *err);

/** Public: interpolate image (as FArray) */
void ObitGPUFInterpolateImage (ObitGPUFInterpolate *in, ObitFArray *inArray,
			       ObitFArray *outArray, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitGPUFInterpolateImageFP) (ObitGPUFInterpolate *in, ObitFArray *inArray,
					    ObitFArray *outArray, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitGPUFInterpolateGetClass (void);

/** Public: Copy (deep) constructor. */
ObitGPUFInterpolate* ObitGPUFInterpolateCopy  (ObitGPUFInterpolate *in, 
					       ObitGPUFInterpolate *out, ObitErr *err);

/** Public: Copy structure. */
void ObitGPUFInterpolateClone (ObitGPUFInterpolate *in, 
			       ObitGPUFInterpolate *out, ObitErr *err);


/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitGPUFInterpolateClassDef.h"
} ObitGPUFInterpolateClassInfo; 

#endif /* OBITFGPUSKYMODEL_H */ 
