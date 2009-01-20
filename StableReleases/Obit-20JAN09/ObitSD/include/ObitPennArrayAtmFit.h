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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITPENNARRAYATMFIT_H 
#define OBITPENNARRAYATMFIT_H 

#include <math.h>
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitOTFArrayGeom.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPennArrayAtmFit.h
 * ObitPennArrayAtmFit Atmospheric model fitting routines for the 
 * Penn Array on the GBT
 * This class is derived from the Obit class.
 *
 */

/*--------------Class definitions-------------------------------------*/
/** ObitPennArrayAtmFit Class structure. */
typedef struct {
#include "ObitPennArrayAtmFitDef.h"   /* this class definition */
} ObitPennArrayAtmFit;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitPennArrayAtmFit
 * returns a ObitPennArrayAtmFit*.
 * in = object to unreference
 */
#define ObitPennArrayAtmFitUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitPennArrayAtmFit.
 * returns a ObitPennArrayAtmFit*.
 * in = object to reference
 */
#define ObitPennArrayAtmFitRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitPennArrayAtmFitIsA(in) ObitIsA (in, ObitPennArrayAtmFitGetClass())

/*---------------Public functions---------------------------*/
/**  Public: Class initializer. */
void ObitPennArrayAtmFitClassInit (void);

/** Public: Constructor. */
ObitPennArrayAtmFit* newObitPennArrayAtmFit (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitPennArrayAtmFitGetClass (void);

/** Public: Constructor from Values. */
ObitPennArrayAtmFit* 
ObitPennArrayAtmFitValue (gchar *name, ObitOTFArrayGeom *geom, olong ncoef, ObitErr *err);

/** Public: Fit model. */
void ObitPennArrayAtmFitFit (ObitPennArrayAtmFit *in, ofloat *data, olong incs,
			     ofloat *coef);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitPennArrayAtmFitClassDef.h"
} ObitPennArrayAtmFitClassInfo; 

#endif /* OBITPENNARRAYATMFIT_H */ 
