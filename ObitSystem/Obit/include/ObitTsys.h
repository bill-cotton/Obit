/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010                                               */
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
#ifndef OBITTSYS_H 
#define OBITTSYS_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitUV.h"
#include "ObitTsys.h"
#include "ObitTableTY.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTsys.h
 *
 * ObitTsys System temperature interpolator class
 *
 * The ObitTsys class provides an interface to a System temperature (TY)
 * table allowing interpolation for a given time and antenna.
 * 
 * \section ObitTsysaccess Creators and Destructors
 * An ObitTsys will usually be created using ObitTsysCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitTsys should always be made using the
 * #ObitTsysRef function which updates the reference count in the object.
 * Then whenever freeing an ObitTsys or changing a pointer, the function
 * #ObitTsysUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitTsys Class structure. */
typedef struct {
#include "ObitTsysDef.h"   /* this class definition */
} ObitTsys;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitTsys
 * returns a ObitTsys*.
 * in = object to unreference
 */
#define ObitTsysUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitTsys.
 * returns a ObitTsys*.
 * in = object to reference
 */
#define ObitTsysRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitTsysIsA(in) ObitIsA (in, ObitTsysGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitTsysClassInit (void);

/** Public: Default Constructor. */
ObitTsys* newObitTsys (gchar* name);

/** Public: Create/initialize ObitTsys structures */
ObitTsys* ObitTsysCreate (gchar* name, ObitTableTY *TYTable, 
			  ObitUV *UVData, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitTsys* (*ObitTsysCreateFP) (gchar* name, ObitTableTY *TYTable, 
				       ObitUV *UVData, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitTsysGetClass (void);

/** Public: Copy (deep) constructor. */
ObitTsys* ObitTsysCopy  (ObitTsys *in, ObitTsys *out, ObitErr *err);

/** Public: Copy structure. */
void ObitTsysClone (ObitTsys *in, ObitTsys *out, ObitErr *err);

/** Public: Interpolate Tsys values. */
void ObitTsysReport (ObitTsys *in, ofloat time, olong ant, olong suba,
		     ofloat *Tsys1, ofloat *Tsys2,
		     ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitTsysClassDef.h"
} ObitTsysClassInfo; 

#endif /* OBITTSYS_H */ 
