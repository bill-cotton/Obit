/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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
#ifndef OBITALMACALATM_H
#define OBITALMACALATM_H

#include "Obit.h"
#include "ObitErr.h"
#include "ObitFile.h"
#include "ObitUVDesc.h"
#include "ObitSDMData.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitALMACalAtm.h
 *
 * This class accesses data in the ALMA binary calAtmosphere format
 *
 * \section ObitALMACalAtmaccess Creators and Destructors
 * An ObitALMACalAtm will usually be created using ObitALMACalAtmCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitALMACalAtm should always be made using the
 * #ObitALMACalAtmRef function which updates the reference count in the object.
 * Then whenever freeing an ObitALMACalAtm or changing a pointer, the function
 * #ObitALMACalAtmUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitESPMIMEType
 * Enum for xml mime type
 */
enum obitACAMIMEType {
  ACAMIMEType_Unknown,
  ACAMIMEType_sdmDataHeader,       /* scan header */
  ACAMIMEType_desc,                /* integration header */
  ACAMIMEType_crossData,           /* cross correlation data */
  ACAMIMEType_autoData,            /* auto correlation data */
  ACAMIMEType_flags,               /* flag data */
  ACAMIMEType_actualTimes,         /* actual time data */
  ACAMIMEType_actualDurations,     /* actual duration data */
  ACAMIMEType_weights,             /* weight data */
  ACAMIMEType_EOF                  /* End of file */
}; /* end enum obitACAMIMEType */
/** typedef for enum for ACAMIMEType. */
typedef enum obitACAMIMEType ObitACAMIMEType;

/**
 * \enum obitALMACalAtmType
 * Enum for data type
 */
enum obitALMACalAtmType {
  ALMACalAtmType_Unknown,
  ALMACalAtmType_INT16_TYPE,
  ALMACalAtmType_INT32_TYPE,
  ALMACalAtmType_INT64_TYPE,
  ALMACalAtmType_FLOAT32_TYPE,
  ALMACalAtmType_FLOAT64_TYPE
}; /* end enum obitALMACalAtmType */
/** typedef for enum for ALMACalAtmType. */
typedef enum obitALMACalAtmType ObitALMACalAtmType;

/**
 * \enum obitACAEndian
 * Enum for endianness of data
 */
enum obitACAEndian {
  ACAEndian_Big,
  ACAEndian_Little
}; /* end enum obitACAEndian */
/** typedef for enum for ACAEndian. */
typedef enum obitACAEndian ObitACAEndian;

/*----------------- Macroes ---------------------------*/
/** Granularity of buffer operations (frame size) */
#define ACABUFFERSIZE 65536
/** Number of frames in buffer */
#define ACABUFFERFRAMES 2
/** 
 * Macro to unreference (and possibly destroy) an ObitALMACalAtm
 * returns a ObitALMACalAtm*.
 * in = object to unreference
 */
#define ObitALMACalAtmUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitALMACalAtm.
 * returns a ObitALMACalAtm*.
 * in = object to reference
 */
#define ObitALMACalAtmRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitALMACalAtmIsA(in) ObitIsA (in, ObitALMACalAtmGetClass())

/*--------------Class definitions-------------------------------------*/
/** ObitALMACalAtm Class structures. */
typedef struct {
#include "ObitALMACalAtmDef.h"   /* this class definition */
} ObitALMACalAtm;

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitALMACalAtmClassInit (void);

/** Public: Default Constructor. */
ObitALMACalAtm* newObitALMACalAtm (gchar* name);

/** Public: Create/initialize ObitALMACalAtm structures */
ObitALMACalAtm* ObitALMACalAtmCreate (gchar* name, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitALMACalAtm* (*ObitALMACalAtmCreateFP) (gchar* name, 
						       ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitALMACalAtmGetClass (void);

/** Public: Copy (deep) constructor. */
ObitALMACalAtm* ObitALMACalAtmCopy  (ObitALMACalAtm *in, ObitALMACalAtm *out, ObitErr *err);

/** Public: Copy structure. */
void ObitALMACalAtmClone (ObitALMACalAtm *in, ObitALMACalAtm *out, ObitErr *err);

/** Public: Initialize file */
void ObitALMACalAtmInitFile (ObitALMACalAtm *in, gchar *DataFile, ObitErr *err);

/** Public: Fill Buffer */
ObitIOCode ObitALMACalAtmFillBuffer (ObitALMACalAtm *in, ObitErr *err);

/** Public: Get number of rows */
olong ObitALMACalAtmGetNrow (ObitALMACalAtm *in);

/** Public: Get visibility record */
olong ObitALMACalAtmGetRow (ObitALMACalAtm *in, ASDMcalAtmosphereRow *row, 
			    ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitALMACalAtmClassDef.h"
} ObitALMACalAtmClassInfo; 

#endif /* OBITALMACALATM_H */ 
