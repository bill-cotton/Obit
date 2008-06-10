/* $Id$  */
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
#ifndef OBITIMAGEFITDATA_H 
#define OBITIMAGEFITDATA_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitImage.h"
#include "ObitFitRegion.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageFitData.h
 * ObitImageFitData Stores fitting data and provides evaluation routines
 * 
 * This class is derived from the #Obit class.
 *
 * \section ObitImageFitDataaccess Creators and Destructors
 * An ObitImageFitData will usually be created using ObitImageFitDataCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitImageFitData should always be made using the
 * #ObitImageFitDataRef function which updates the reference count in the object.
 * Then whenever freeing an ObitImageFitData or changing a pointer, the function
 * #ObitImageFitDataUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** function pointer for evaluation routine used by solver */
typedef void (*ObitImageFitDataFuncFP) 
     (gpointer data, odouble* p, odouble *f, odouble* grad, 
      olong iflag);

/** ObitImageFitData Class structure. */
typedef struct {
#include "ObitImageFitDataDef.h"   /* this class definition */
} ObitImageFitData;

/*-------------- enumerations -------------------------------------*/
/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitImageFitData
 * returns a ObitImageFitData*.
 * in = object to unreference
 */
#define ObitImageFitDataUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitImageFitData.
 * returns a ObitImageFitData*.
 * in = object to reference
 */
#define ObitImageFitDataRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitImageFitDataIsA(in) ObitIsA (in, ObitImageFitDataGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitImageFitDataClassInit (void);

/** Public: Default Constructor. */
ObitImageFitData* newObitImageFitData (gchar* name);

/** Public: Create/initialize ObitImageFitData structures */
ObitImageFitData* 
ObitImageFitDataCreate (gchar* name, ObitFitRegion *reg, ObitInfoList *bounds, 
			ObitImage *image, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitImageFitData* 
(*ObitImageFitDataCreateFP) (gchar* name, ObitFitRegion *reg, 
			     ObitInfoList *bounds, ObitImage *image, 
			     ObitErr *err);

/** Public: Copy contents to ObitFitRegion */
void ObitImageFitData2Reg (ObitImageFitData* data, ObitFitRegion *reg);
typedef void (*ObitImageFitData2RegFP) 
     (ObitImageFitData* data, ObitFitRegion *reg);

/** Public: Determine errors for Gaussian model */
void ObitImageFitDataGaussErr (ofloat peak, ofloat major, ofloat minor, ofloat posang, 
			       ofloat irms, ofloat* beam, 
			       ofloat* epeak, ofloat* errra, ofloat* errdec, 
			       ofloat* errmaj, ofloat* errmin, ofloat* errpa);

/** Public: ClassInfo pointer */
gconstpointer ObitImageFitDataGetClass (void);
/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitImageFitDataClassDef.h"
} ObitImageFitDataClassInfo; 

#endif /* OBITIMAGEFITDATA_H */ 
