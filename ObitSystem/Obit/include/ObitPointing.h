/* $Id: ObitPointing.h 467 2013-12-20 14:10:44Z bill.cotton $  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2016                                               */
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
#ifndef OBITPOINTING_H
#define OBITPOINTING_H

#include "Obit.h"
#include "ObitErr.h"
#include "ObitFile.h"
#include "ObitUVDesc.h"
#include "ObitSDMData.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPointing.h
 *
 * This class accesses data in ASDM Pointing table binary format
 *
 * \section ObitPointingaccess Creators and Destructors
 * An ObitPointing will usually be created using ObitPointingCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitPointing should always be made using the
 * #ObitPointingRef function which updates the reference count in the object.
 * Then whenever freeing an ObitPointing or changing a pointer, the function
 * #ObitPointingUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitPointingType
 * Enum for data type
 */
enum obitPointingType {
  PointingType_Unknown,
  PointingType_INT16_TYPE,
  PointingType_INT32_TYPE,
  PointingType_INT64_TYPE,
  PointingType_FLOAT32_TYPE,
  PointingType_FLOAT64_TYPE
}; /* end enum obitPointingType */
/** typedef for enum for PointingType. */
typedef enum obitPointingType ObitPointingType;

/**
 * \enum obitPointingEndian
 * Enum for endianness of data
 */
enum obitPointingEndian {
  PointingEndian_Big,
  PointingEndian_Little
}; /* end enum obitPointingEndian */
/** typedef for enum for PointinEndian. */
typedef enum obitPointingEndian ObitPointingEndian;

/*----------------- Macroes ---------------------------*/
/** Granularity of buffer operations (frame size) */
#define POINTINGBUFFERSIZE 32768*8
/** Number of frames in buffer */
#define POINTINGBUFFERFRAMES 2
/** 
 * Macro to unreference (and possibly destroy) an ObitPointing
 * returns a ObitPointing*.
 * in = object to unreference
 */
#define ObitPointingUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitPointing.
 * returns a ObitPointing*.
 * in = object to reference
 */
#define ObitPointingRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitPointingIsA(in) ObitIsA (in, ObitPointingGetClass())

/*--------------Class definitions-------------------------------------*/
/** ObitPointing Class structures. */
typedef struct {
#include "ObitPointingDef.h"   /* this class definition */
} ObitPointing;

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitPointingClassInit (void);

/** Public: Default Constructor. */
ObitPointing* newObitPointing (gchar* name);

/** Public: Create/initialize ObitPointing structures */
ObitPointing* ObitPointingCreate (gchar* name, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitPointing* (*ObitPointingCreateFP) (gchar* name, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitPointingGetClass (void);

/** Public: Copy (deep) constructor. */
ObitPointing* ObitPointingCopy  (ObitPointing *in, ObitPointing *out, ObitErr *err);

/** Public: Copy structure. */
void ObitPointingClone (ObitPointing *in, ObitPointing *out, ObitErr *err);

/** Public: Initialize file */
void ObitPointingInitFile (ObitPointing *in, gchar *DataFile, ObitErr *err);

/** Public: Fill Buffer */
ObitIOCode ObitPointingFillBuffer (ObitPointing *in, ObitErr *err);

/** Public: Get number of rows */
olong ObitPointingGetNrow (ObitPointing *in);

/** Public: Get pointing record */
olong ObitPointingGetRow (ObitPointing *in, ASDMPointingRow *row, 
			  ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
  #include "ObitPointingClassDef.h"
} ObitPointingClassInfo; 

#endif /* OBITPOINTING_H */ 
