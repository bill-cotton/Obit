/* $Id$        */
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
#include "ObitCUDASkyModel.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGPUSkyModel.h
 *
 * ObitGPUSkyModel GPU enhanced base sky model class
 * Uses functions in ObitCUDASkyModel
 *
 * Sky models are calculated using a GPU.
 * 
 * \section ObitGPUSkyModelaccess Creators and Destructors
 * An ObitGPUSkyModel will usually be created using ObitGPUSkyModelCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitGPUSkyModel should always be made using the
 * #ObitGPUSkyModelRef function which updates the reference count in the object.
 * Then whenever freeing an ObitGPUSkyModel or changing a pointer, the function
 * #ObitGPUSkyModelUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitGPUSkyModel Class structure. */
typedef struct {
#include "ObitGPUSkyModelDef.h"   /* this class definition */
} ObitGPUSkyModel;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitGPUSkyModel
 * returns a ObitGPUSkyModel*.
 * in = object to unreference
 */
#define ObitGPUSkyModelUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitGPUSkyModel.
 * returns a ObitGPUSkyModel*.
 * in = object to reference
 */
#define ObitGPUSkyModelRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGPUSkyModelIsA(in) ObitIsA (in, ObitGPUSkyModelGetClass())

/*---------------------------------- Structures ----------------------------*/
/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitGPUSkyModelClassInit (void);

/** Public: Default Constructor. */
ObitGPUSkyModel* newObitGPUSkyModel (gchar* name);

/** Public: Create/initialize ObitGPUSkyModel structures */
ObitGPUSkyModel* ObitGPUSkyModelCreate (gchar* name, gchar *type);
/** Typedef for definition of class pointer structure */
typedef ObitGPUSkyModel* (*ObitGPUSkyModelCreateFP) (gchar* name, gchar *type);

/** Public: ClassInfo pointer */
gconstpointer ObitGPUSkyModelGetClass (void);

/** Public: Copy (deep) constructor. */
ObitGPUSkyModel* ObitGPUSkyModelCopy  (ObitGPUSkyModel *in, ObitGPUSkyModel *out, ObitErr *err);

/** Public: Copy structure. */
void ObitGPUSkyModelClone (ObitGPUSkyModel *in, ObitGPUSkyModel *out, ObitErr *err);

/* Public: Initialize DFT Model */
void ObitGPUSkyModelDFTInit (ObitGPUSkyModel *in, Obit *skyModel,
			     ObitUV *uvdata, ObitErr *err);
typedef void (*ObitGPUSkyModelDFTInitFP) (ObitGPUSkyModel *in, Obit *skyModel, 
					  ObitUV *uvdata, ObitErr *err);

/* Public: Set DFT sky model Model */
void ObitGPUSkyModelDFTSetMod (ObitGPUSkyModel *in, Obit *skyModel,
			       ObitFArray *model, ObitErr *err);
typedef void (*ObitGPUSkyModelDFTSetModFP) (ObitGPUSkyModel *in, Obit *skyModel,
					    ObitFArray *model, ObitErr *err);

/* Public: Calculate DFT Model */
void ObitGPUSkyModelDFTCalc (ObitGPUSkyModel *in, ObitUV *uvdata, ObitErr *err);
typedef void (*ObitGPUSkyModelDFTCalcFP) (ObitGPUSkyModel *in, ObitUV *uvdata, ObitErr *err);

/* Public: Shutdown DFT Model */
void ObitGPUSkyModelDFTShutdown (ObitGPUSkyModel *in, ObitUV *uvdata, ObitErr *err);
typedef void (*ObitGPUSkyModelDFTShutdownFP) (ObitGPUSkyModel *in, ObitUV *uvdata, ObitErr *err);


/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitGPUSkyModelClassDef.h"
} ObitGPUSkyModelClassInfo; 

#endif /* OBITFGPUSKYMODEL_H */ 
