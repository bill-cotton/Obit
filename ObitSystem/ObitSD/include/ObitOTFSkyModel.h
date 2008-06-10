/* $Id: ObitOTFSkyModel.h,v 1.4 2005/10/06 19:33:28 bcotton Exp $  */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITOTFSKYMODEL_H 
#define OBITOTFSKYMODEL_H 

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitTableSkyModel.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFSkyModel.h
 * OTF Sky model
 *
 * This class is for creating and manipulating sky model objects
 * for OTF data.
 * 
 * \section ObitOTFSkyModelaccess Creators and Destructors
 * An ObitOTFSkyModel will usually be created using ObitOTFSkyModelCreate which allows 
 * specifying a name for the object as well as dimensionality of the array.
 *
 * A copy of a pointer to an ObitOTFSkyModel should always be made using the
 * #ObitOTFSkyModelRef function which updates the reference count in the object.
 * Then whenever freeing an ObitOTFSkyModel or changing a pointer, the function
 * #ObitOTFSkyModelUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitOTFProj
 * enum for OTF projection types.
 */
enum obitOTFProj {
  /** -SIN  projection */
  OBIT_OTF_SIN = 0, 
  /** -ARC projection */
  OBIT_OTF_ARC, 
  /** -TAN  projection */
  OBIT_OTF_TAN 
}; /* end enum obitIOType */

/** typedef for enum for OTF projection types. */
typedef enum obitOTFProj ObitOTFProj;

/*--------------Class definitions-------------------------------------*/
/** ObitOTFSkyModel Class structure. */
typedef struct {
#include "ObitOTFSkyModelDef.h"   /* this class definition */
} ObitOTFSkyModel;

/*----------------------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitOTFSkyModel
 * returns a ObitOTFSkyModel*.
 * in = object to unreference
 */
#define ObitOTFSkyModelUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitOTFSkyModel.
 * returns a ObitOTFSkyModel*.
 * in = object to reference
 */
#define ObitOTFSkyModelRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitOTFSkyModelIsA(in) ObitIsA (in, ObitOTFSkyModelGetClass())

/*---------------Public functions---------------------------*/
/** Public : Class initializer. */
void ObitOTFSkyModelClassInit (void);

/** Public: Default Constructor. */
ObitOTFSkyModel* newObitOTFSkyModel (gchar* name);

/** Public: Create/initialize ObitOTFSkyModel structures */
ObitOTFSkyModel* ObitOTFSkyModelCreate (olong ndetect);
/** Typedef for definition of class pointer structure */
typedef void (*ObitOTFSkyModelCreateFP) (olong ndetect);

/** Public: ClassInfo pointer */
gconstpointer ObitOTFSkyModelGetClass (void);

/** Public: Copy (deep) constructor. */
ObitOTFSkyModel* 
ObitOTFSkyModelCopy  (ObitOTFSkyModel *in, ObitOTFSkyModel *out, ObitErr *err);

/** Public:  Read Table from disk */
ObitIOCode ObitOTFSkyModelRead (ObitOTFSkyModel **in, ObitTableSkyModel *table, ObitErr *err);

/** Public:  Write Table to disk */
ObitIOCode ObitOTFSkyModelWrite (ObitOTFSkyModel *in, ObitTableSkyModel *table, ObitErr *err);

/** Public: Determine Projection type from a string */
ObitOTFProj ObitOTFSkyModelProj (gchar *string);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitOTFSkyModelClassDef.h"
} ObitOTFSkyModelClassInfo; 

#endif /* OBITOTFSKYMODEL_H */ 
