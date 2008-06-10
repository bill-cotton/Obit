/* $Id: ObitFitModel.h,v 1.2 2007/08/31 17:24:48 bcotton Exp $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2008                                          */
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
#ifndef OBITFITMODEL_H 
#define OBITFITMODEL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitImage.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFitModel.h
 * ObitFitModel Stores model fit parameters
 *
 * This class is derived from the #Obit class.
 * 
 * \section ObitFitModelaccess Creators and Destructors
 * An ObitFitModel will usually be created using ObitFitModelCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitFitModel should always be made using the
 * #ObitFitModelRef function which updates the reference count in the object.
 * Then whenever freeing an ObitFitModel or changing a pointer, the function
 * #ObitFitModelUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitFitModelCompType
 * enum for Model Component model type
 */
enum obitFitModelCompType {
  /** Point */
  OBIT_FitModel_PointMod,
  /** Gaussian on sky */
  OBIT_FitModel_GaussMod,  
  /** Convolved Gaussian */
  OBIT_FitModel_CGaussMod,  
  /** Uniform sphere */
  OBIT_FitModel_USphereMod, 
  /** Background surface */
  OBIT_FitModel_Background, 
  /** Unknown */
  OBIT_FitModel_Unknown 
}; /* end enum obitFitModelCompType */
/** typedef for enum for ObitFitModelCompType. */
typedef enum obitFitModelCompType ObitFitModelCompType;

/*--------------Class definitions-------------------------------------*/
/** ObitFitModel Class structure. */
typedef struct {
#include "ObitFitModelDef.h"   /* this class definition */
} ObitFitModel;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitFitModel
 * returns a ObitFitModel*.
 * in = object to unreference
 */
#define ObitFitModelUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitFitModel.
 * returns a ObitFitModel*.
 * in = object to reference
 */
#define ObitFitModelRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitFitModelIsA(in) ObitIsA (in, ObitFitModelGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitFitModelClassInit (void);

/** Public: Default Constructor. */
ObitFitModel* newObitFitModel (gchar* name);

/** Public: Create/initialize ObitFitModel structures */
ObitFitModel* 
ObitFitModelCreate (gchar* name, ObitFitModelCompType type, 
		    ofloat Flux, ofloat DeltaX, ofloat DeltaY, 
		    olong nparm, ofloat *parms);
/** Typedef for definition of class pointer structure */
typedef ObitFitModel* 
(*ObitFitModelCreateFP) (gchar* name, ObitFitModelCompType type, 
			 ofloat Flux, ofloat DeltaX, ofloat DeltaY, 
			 olong nparm, ofloat *parms);

/** Public: ClassInfo pointer */
gconstpointer ObitFitModelGetClass (void);

/** Public: Copy (deep) constructor. */
ObitFitModel* ObitFitModelCopy  (ObitFitModel *in, ObitFitModel *out, 
				 ObitErr *err);

/** Public: Copy structure. */
void ObitFitModelClone (ObitFitModel *in, ObitFitModel *out, ObitErr *err);

/** Public: Deconvolve Gaussians. */
olong ObitFitModelDeconGau (ofloat bmaj, ofloat bmin, ofloat bpa, 
			   ofloat ebmaj, ofloat ebmin, ofloat ebpa,
			   ofloat cbmaj, ofloat cbmin, ofloat cbpa,
			   ofloat dgau[3][3]);
/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitFitModelClassDef.h"
} ObitFitModelClassInfo; 

#endif /* OBITFITMODEL_H */ 
