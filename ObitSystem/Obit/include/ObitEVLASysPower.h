/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2011                                               */
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
#ifndef OBITEVLASYSPOWER_H 
#define OBITEVLASYSPOWER_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitFile.h"
#include "ObitUVDesc.h"
#include "ObitSDMData.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitEVLASysPower.h
 *
 * This class accesses data in the EVLA binary SysPower format
 *
 * \section ObitEVLASysPoweraccess Creators and Destructors
 * An ObitEVLASysPower will usually be created using ObitEVLASysPowerCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitEVLASysPower should always be made using the
 * #ObitEVLASysPowerRef function which updates the reference count in the object.
 * Then whenever freeing an ObitEVLASysPower or changing a pointer, the function
 * #ObitEVLASysPowerUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitESPMIMEType
 * Enum for xml mime type
 */
enum obitESPMIMEType {
  ESPMIMEType_Unknown,
  ESPMIMEType_sdmDataHeader,       /* scan header */
  ESPMIMEType_desc,                /* integration header */
  ESPMIMEType_crossData,           /* cross correlation data */
  ESPMIMEType_autoData,            /* auto correlation data */
  ESPMIMEType_flags,               /* flag data */
  ESPMIMEType_actualTimes,         /* actual time data */
  ESPMIMEType_actualDurations,     /* actual duration data */
  ESPMIMEType_weights,             /* weight data */
  ESPMIMEType_EOF                  /* End of file */
}; /* end enum obitESPMIMEType */
/** typedef for enum for ESPMIMEType. */
typedef enum obitESPMIMEType ObitESPMIMEType;

/**
 * \enum obitEVLASysPowerType
 * Enum for data type
 */
enum obitEVLASysPowerType {
  EVLASysPowerType_Unknown,
  EVLASysPowerType_INT16_TYPE,
  EVLASysPowerType_INT32_TYPE,
  EVLASysPowerType_INT64_TYPE,
  EVLASysPowerType_FLOAT32_TYPE,
  EVLASysPowerType_FLOAT64_TYPE
}; /* end enum obitEVLASysPowerType */
/** typedef for enum for EVLASysPowerType. */
typedef enum obitEVLASysPowerType ObitEVLASysPowerType;

/**
 * \enum obitESPEndian
 * Enum for endianness of data
 */
enum obitESPEndian {
  ESPEndian_Big,
  ESPEndian_Little
}; /* end enum obitESPEndian */
/** typedef for enum for ESPEndian. */
typedef enum obitESPEndian ObitESPEndian;

/*----------------- Macroes ---------------------------*/
/** Granularity of buffer operations (frame size) */
#define ESPBUFFERSIZE 2048
/** Number of frames in buffer */
#define ESPBUFFERFRAMES 2
/** 
 * Macro to unreference (and possibly destroy) an ObitEVLASysPower
 * returns a ObitEVLASysPower*.
 * in = object to unreference
 */
#define ObitEVLASysPowerUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitEVLASysPower.
 * returns a ObitEVLASysPower*.
 * in = object to reference
 */
#define ObitEVLASysPowerRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitEVLASysPowerIsA(in) ObitIsA (in, ObitEVLASysPowerGetClass())

/*--------------Class definitions-------------------------------------*/
/** ObitEVLASysPower Class structures. */
typedef struct {
#include "ObitEVLASysPowerDef.h"   /* this class definition */
} ObitEVLASysPower;

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitEVLASysPowerClassInit (void);

/** Public: Default Constructor. */
ObitEVLASysPower* newObitEVLASysPower (gchar* name);

/** Public: Create/initialize ObitEVLASysPower structures */
ObitEVLASysPower* ObitEVLASysPowerCreate (gchar* name, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitEVLASysPower* (*ObitEVLASysPowerCreateFP) (gchar* name, 
						       ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitEVLASysPowerGetClass (void);

/** Public: Copy (deep) constructor. */
ObitEVLASysPower* ObitEVLASysPowerCopy  (ObitEVLASysPower *in, ObitEVLASysPower *out, ObitErr *err);

/** Public: Copy structure. */
void ObitEVLASysPowerClone (ObitEVLASysPower *in, ObitEVLASysPower *out, ObitErr *err);

/** Public: Initialize file */
void ObitEVLASysPowerInitFile (ObitEVLASysPower *in, gchar *DataFile, ObitErr *err);

/** Public: Fill Buffer */
ObitIOCode ObitEVLASysPowerFillBuffer (ObitEVLASysPower *in, ObitErr *err);

/** Public: Get number of rows */
olong ObitEVLASysPowerGetNrow (ObitEVLASysPower *in);

/** Public: Get visibility record */
olong ObitEVLASysPowerGetRow (ObitEVLASysPower *in, ASDMSysPowerRow *row, 
			      ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitEVLASysPowerClassDef.h"
} ObitEVLASysPowerClassInfo; 

#endif /* OBITFEVLASYSPOWER_H */ 
