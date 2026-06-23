/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2026                                               */
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
#ifndef OBITGPUBEAMMODEL_H 
#define OBITGPUBEAMMODEL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitUV.h"
#include "ObitGPUBeamInterp.h"
#include "ObitCUDABeamModel.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGPUBeamModel.h
 *
 * ObitGPUBeamModel GPU enhanced base sky model class with Beam corrections
 * Uses functions in ObitCUDASkyModel
 *
 * Sky models are calculated using a GPU.
 * 
 * \section ObitGPUBeamModelaccess Creators and Destructors
 * An ObitGPUBeamModel will usually be created using ObitGPUBeamModelCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitGPUBeamModel should always be made using the
 * #ObitGPUBeamModelRef function which updates the reference count in the object.
 * Then whenever freeing an ObitGPUBeamModel or changing a pointer, the function
 * #ObitGPUBeamModelUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitGPUBeamModel Class structure. */
typedef struct {
#include "ObitGPUBeamModelDef.h"   /* this class definition */
} ObitGPUBeamModel;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitGPUBeamModel
 * returns a ObitGPUBeamModel*.
 * in = object to unreference
 */
#define ObitGPUBeamModelUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitGPUBeamModel.
 * returns a ObitGPUBeamModel*.
 * in = object to reference
 */
#define ObitGPUBeamModelRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGPUBeamModelIsA(in) ObitIsA (in, ObitGPUBeamModelGetClass())

/*---------------------------------- Structures ----------------------------*/
/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitGPUBeamModelClassInit (void);

/** Public: Default Constructor. */
ObitGPUBeamModel* newObitGPUBeamModel (gchar* name);

/** Public: Create/initialize ObitGPUBeamModel structures */
ObitGPUBeamModel* ObitGPUBeamModelCreate (gchar* name, gchar *type);
/** Typedef for definition of class pointer structure */
typedef ObitGPUBeamModel* (*ObitGPUBeamModelCreateFP) (gchar* name, gchar *type);

/** Public: ClassInfo pointer */
gconstpointer ObitGPUBeamModelGetClass (void);

/** Public: Copy (deep) constructor. */
ObitGPUBeamModel* ObitGPUBeamModelCopy  (ObitGPUBeamModel *in, ObitGPUBeamModel *out, ObitErr *err);

/** Public: Copy structure. */
void ObitGPUBeamModelClone (ObitGPUBeamModel *in, ObitGPUBeamModel *out, ObitErr *err);

/* Public: Initialize DFT Model */
void ObitGPUBeamModelDFTInit (ObitGPUBeamModel *in, Obit *skyModel,
			     ObitUV *uvdata, ofloat delta_PA,  ObitErr *err);
typedef void (*ObitGPUBeamModelDFTInitFP) (ObitGPUBeamModel *in, Obit *skyModel, 
					   ObitUV *uvdata, ofloat delta_PA,
					   ObitErr *err);

/* Public: Set DFT sky model Model */
void ObitGPUBeamModelDFTSetMod (ObitGPUBeamModel *in, Obit *skyModel,
			       ObitFArray *model, ObitErr *err);
typedef void (*ObitGPUBeamModelDFTSetModFP) (ObitGPUBeamModel *in, Obit *skyModel,
					    ObitFArray *model, ObitErr *err);

/* public: setup antenna beams */
void ObitGPUBeamModelDFTSetBeam (ObitGPUBeamModel *in, Obit *SkyModel, ObitErr *err);
typedef void (*ObitGPUBeamModelDFTSetBeamFP) (ObitGPUBeamModel *in, Obit *skyModel,
					      ObitErr *err);


/* Public: Calculate DFT Model */
void ObitGPUBeamModelDFTCalc (ObitGPUBeamModel *in, ObitUV *uvdata, ObitErr *err);
typedef void (*ObitGPUBeamModelDFTCalcFP) (ObitGPUBeamModel *in, ObitUV *uvdata, ObitErr *err);

/* Public: Shutdown DFT Model */
void ObitGPUBeamModelDFTShutdown (ObitGPUBeamModel *in, ObitUV *uvdata, ObitErr *err);
typedef void (*ObitGPUBeamModelDFTShutdownFP) (ObitGPUBeamModel *in, ObitUV *uvdata, ObitErr *err);


/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitGPUBeamModelClassDef.h"
} ObitGPUBeamModelClassInfo; 

#endif /* OBITFGPUBEAMMODEL_H */ 
