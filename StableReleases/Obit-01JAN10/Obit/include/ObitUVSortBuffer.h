/* $Id:  $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009                                               */
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
#ifndef OBITUVSORTBUFFER_H 
#define OBITUVSORTBUFFER_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitUV.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/

/**
 * \file ObitUVSortBuffer.h
 *
 * ObitUVSortBuffer Sorting buffer for UV data
 *
 * The ObitUVSortBuffer assists in the sorting of UV data into time order 
 * by providing an  ObitUV data buffer with sorting and other 
 * manipulation facilities. 
 * 
 * \section ObitUVSortBufferaccess Creators and Destructors
 * An ObitUVSortBuffer will usually be created using ObitUVSortBufferCreate 
 * which allows specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitUVSortBuffer should always be made using the
 * #ObitUVSortBufferRef function which updates the reference count in the object.
 * Then whenever freeing an ObitUVSortBuffer or changing a pointer, the function
 * #ObitUVSortBufferUnref will decrement the reference count and destroy the 
 * object when the reference count hits 0.
 * There is no explicit destructor.
 */
/*-------------------- structure -------------------------------------*/
/** Equivalence for Sort key arrays */
  union ObitUVSortEquiv { 
    olong   itg;
    ofloat  flt[2];
  };
/** Sort Index plus sort key */
typedef struct {
  union ObitUVSortEquiv index;
  ofloat key[2];
} ObitUVSortStruct;

/*--------------Class definitions-------------------------------------*/
/** ObitUVSortBuffer Class structure. */
typedef struct {
#include "ObitUVSortBufferDef.h"   /* this class definition */
} ObitUVSortBuffer;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitUVSortBuffer
 * returns a ObitUVSortBuffer*.
 * in = object to unreference
 */
#define ObitUVSortBufferUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitUVSortBuffer.
 * returns a ObitUVSortBuffer*.
 * in = object to reference
 */
#define ObitUVSortBufferRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitUVSortBufferIsA(in) ObitIsA (in, ObitUVSortBufferGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitUVSortBufferClassInit (void);

/** Public: Default Constructor. */
ObitUVSortBuffer* newObitUVSortBuffer (gchar* name);

/** Public: Create/initialize ObitUVSortBuffer structures */
ObitUVSortBuffer* ObitUVSortBufferCreate (gchar* name, ObitUV *inUV, 
					  olong nvis, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitUVSortBuffer* (*ObitUVSortBufferCreateFP) (gchar* name, 
						       ObitUV *inUV, 
						       olong nvis, 
						       ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitUVSortBufferGetClass (void);

/** Public: Copy (deep) constructor. */
ObitUVSortBuffer* ObitUVSortBufferCopy  (ObitUVSortBuffer *in, 
					 ObitUVSortBuffer *out, 
					 ObitErr *err);

/** Public: Copy structure. */
void ObitUVSortBufferClone (ObitUVSortBuffer *in, ObitUVSortBuffer *out, 
			    ObitErr *err);


/** Public: Add a visibility. */
void ObitUVSortBufferAddVis (ObitUVSortBuffer *in,  ofloat *vis, 
			     ofloat lastTime, ObitErr *err);
typedef void (*ObitUVSortBufferAddVisFP) (ObitUVSortBuffer *in,  ofloat *vis, 
					  ofloat lastTime, ObitErr *err);

/** Public: Sort buffer. */
void ObitUVSortBufferSort (ObitUVSortBuffer *in,  ObitErr *err);
typedef void (*ObitUVSortBufferSortFP) (ObitUVSortBuffer *in,  
					  ObitErr *err);

/** Public: Flush buffer. */
void ObitUVSortBufferFlush (ObitUVSortBuffer *in,  ObitErr *err);
typedef void (*ObitUVSortBufferFlushFP) (ObitUVSortBuffer *in,  
					  ObitErr *err);
/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitUVSortBufferClassDef.h"
} ObitUVSortBufferClassInfo; 

#endif /* OBITFUVSORTBUFFER_H */ 
