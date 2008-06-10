/* $Id$ */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITFINTERPOLATE_H 
#define OBITFINTERPOLATE_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitImageDesc.h"
#include "ObitFArray.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFInterpolate.h
 * ObitFInterpolate does Lagrangian interpolation of positions in an ObitFArray
 *
 * This class is derived from the #Obit class.
 *
 */

/*--------------Class definitions-------------------------------------*/
/** ObitFInterpolate Class structure. */
typedef struct {
#include "ObitFInterpolateDef.h"   /* this class definition */
} ObitFInterpolate;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitFInterpolate
 * returns a ObitFInterpolate*.
 * in = object to unreference
 */
#define ObitFInterpolateUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitFInterpolate.
 * returns a ObitFInterpolate*.
 * in = object to reference
 */
#define ObitFInterpolateRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitFInterpolateIsA(in) ObitIsA (in, ObitFInterpolateGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitFInterpolateClassInit (void);

/** Public: Constructor. */
ObitFInterpolate* newObitFInterpolate (gchar* name);

/** Public: Constructor from value. */
ObitFInterpolate* 
newObitFInterpolateCreate (gchar* name, ObitFArray *array, ObitImageDesc *desc, 
			   olong hwidth);

/** Public: ClassInfo pointer */
gconstpointer ObitFInterpolateGetClass (void);

/** Public: Copy (deep) constructor. */
ObitFInterpolate* ObitFInterpolateCopy  (ObitFInterpolate *in, ObitFInterpolate *out, 
			   ObitErr *err);

/** Public: Copy (shallow) constructor. */
ObitFInterpolate* ObitFInterpolateClone (ObitFInterpolate *in, ObitFInterpolate *out);

/** Public: Replace member ObitFArray*/
void ObitFInterpolateReplace (ObitFInterpolate *in, ObitFArray *newArray);

/** Public: Interpolate Pixel in 2D array */
ofloat ObitFInterpolatePixel (ObitFInterpolate *in, ofloat *pixel, ObitErr *err);

/** Public: Interpolate value in 1- array */
ofloat ObitFInterpolate1D (ObitFInterpolate *in, ofloat pixel);

/** Public: Interpolate Position in 2D array */
ofloat ObitFInterpolatePosition (ObitFInterpolate *in, odouble *coord, ObitErr *err);

/** Public: Interpolate Offset in 2D array */
ofloat ObitFInterpolateOffset (ObitFInterpolate *in, odouble *offset, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitFInterpolateClassDef.h"
} ObitFInterpolateClassInfo; 

#endif /* OBITFINTERPOLATE_H */ 
