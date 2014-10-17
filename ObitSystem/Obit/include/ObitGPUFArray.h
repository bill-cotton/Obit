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
#ifndef OBITGPUFARRAY_H 
#define OBITGPUFARRAY_H 

#include "Obit.h"
#include "ObitErr.h"
#include "CUDAFArray.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGPUFArray.h
 *
 * ObitGPUFArray GPU enhanced FArray interpolation class
 * Uses functions in ObitCUDAFArray
 *
 * ObitGPUFArrays consist of an array in locked host memory and corresponding
 * allocations in the GPU.  Data movement functions are provided.
 * Communications to the GPU is via CUDAFArray.
 * 
 * \section ObitGPUFArrayaccess Creators and Destructors
 * An ObitGPUFArray will usually be created using ObitGPUFArrayCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitGPUFArray should always be made using the
 * #ObitGPUFArrayRef function which updates the reference count in the object.
 * Then whenever freeing an ObitGPUFArray or changing a pointer, the function
 * #ObitGPUFArrayUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitGPUFArray Class structure. */
typedef struct {
#include "ObitGPUFArrayDef.h"   /* this class definition */
} ObitGPUFArray;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitGPUFArray
 * returns a ObitGPUFArray*.
 * in = object to unreference
 */
#define ObitGPUFArrayUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitGPUFArray.
 * returns a ObitGPUFArray*.
 * in = object to reference
 */
#define ObitGPUFArrayRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGPUFArrayIsA(in) ObitIsA (in, ObitGPUFArrayGetClass())

/*---------------------------------- Structures ----------------------------*/
/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitGPUFArrayClassInit (void);

/** Public: Default Constructor. */
ObitGPUFArray* newObitGPUFArray (gchar* name);

/** Public: Create/initialize ObitGPUFArray structures */
ObitGPUFArray* ObitGPUFArrayCreate (gchar* name, olong ndim, olong *naxis);
/** Typedef for definition of class pointer structure */
typedef ObitGPUFArray* (*ObitGPUFArrayCreateFP) (gchar* name, 
						 olong ndim, olong *naxis);

/** Public: ClassInfo pointer */
gconstpointer ObitGPUFArrayGetClass (void);

/** Public: Copy (deep) constructor. */
ObitGPUFArray* ObitGPUFArrayCopy  (ObitGPUFArray *in, ObitGPUFArray *out, ObitErr *err);

/** Public: Copy structure. */
void ObitGPUFArrayClone (ObitGPUFArray *in, ObitGPUFArray *out, ObitErr *err);

/* Public: Copy host array to data */
void ObitGPUFArrayToData (ObitGPUFArray *in, ofloat* data, ObitErr *err);
typedef void (*ObitGPUFArrayToDataFP) (ObitGPUFArray *in, ofloat* data, 
				       ObitErr *err);

/* Public:  Copy from data to host array */
void ObitGPUFArrayFromData (ObitGPUFArray *in, ofloat* data, ObitErr *err);
typedef void (*ObitGPUFArrayFromDataFP) (ObitGPUFArray *in, ofloat* data, 
					 ObitErr *err);

/* Public:  Copy array data to GPU */
void ObitGPUFArrayToGPU (ObitGPUFArray *in, ObitErr *err);
typedef void (*ObitGPUFArrayToGPUFP) (ObitGPUFArray *in, ObitErr *err);

/* Public: Copy array data to Host */
void ObitGPUFArrayToHost (ObitGPUFArray *in, ObitErr *err);
typedef void (*ObitGPUFArrayToHostFP) (ObitGPUFArray *in, ObitErr *err);


/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitGPUFArrayClassDef.h"
} ObitGPUFArrayClassInfo; 

#endif /* OBITFGPUFARRAY_H */ 
