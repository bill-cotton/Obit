/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                          */
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
#ifndef OBITCINTERPOLATE_H 
#define OBITCINTERPOLATE_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitImageDesc.h"
#include "ObitCArray.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitCInterpolate.h
 * ObitCInterpolate does Lagrangian interpolation of positions in an ObitCArray
 *
 * This class is derived from the #Obit class.
 *
 */

/*--------------Class definitions-------------------------------------*/
/** ObitCInterpolate Class structure. */
typedef struct {
#include "ObitCInterpolateDef.h"   /* this class definition */
} ObitCInterpolate;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitCInterpolate
 * returns a ObitCInterpolate*.
 * in = object to unreference
 */
#define ObitCInterpolateUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitCInterpolate.
 * returns a ObitCInterpolate*.
 * in = object to reference
 */
#define ObitCInterpolateRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitCInterpolateIsA(in) ObitIsA (in, ObitCInterpolateGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitCInterpolateClassInit (void);

/** Public: Constructor. */
ObitCInterpolate* newObitCInterpolate (gchar* name);

/** Public: Constructor from value. */
ObitCInterpolate* 
newObitCInterpolateCreate (gchar* name, ObitCArray *array, ObitImageDesc *desc, 
			   ofloat OSX, ofloat OSY, olong numConjCol, olong hwidth, 
			   ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitCInterpolateGetClass (void);

/** Public: Copy (deep) constructor. */
ObitCInterpolate* ObitCInterpolateCopy  (ObitCInterpolate *in, ObitCInterpolate *out, 
			   ObitErr *err);

/** Public: Copy (shallow) constructor. */
ObitCInterpolate* ObitCInterpolateClone (ObitCInterpolate *in, ObitCInterpolate *out);

/** Public: Replace member ObitCArray*/
void ObitCInterpolateReplace (ObitCInterpolate *in, ObitCArray *newArray);

/** Public: Interpolate Pixel in 2D array */
void ObitCInterpolatePixel (ObitCInterpolate *in, ofloat *pixel, ofloat out[2], 
			    ObitErr *err);

/** Public: Interpolate value in 1- array */
void ObitCInterpolate1D (ObitCInterpolate *in, ofloat pixel, ofloat out[2]);

/** Public: Interpolate Position in 2D array */
void ObitCInterpolatePosition (ObitCInterpolate *in, odouble *coord, ofloat out[2], 
			       ObitErr *err);

/** Public: Interpolate Offset in 2D array */
void ObitCInterpolateOffset (ObitCInterpolate *in, ofloat *offset, ofloat out[2],
			     ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitCInterpolateClassDef.h"
} ObitCInterpolateClassInfo; 

#endif /* OBITCINTERPOLATE_H */ 
