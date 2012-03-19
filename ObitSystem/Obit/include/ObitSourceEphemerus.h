/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012                                               */
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
#ifndef OBITSOURCEEPHEMERUS_H 
#define OBITSOURCEEPHEMERUS_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitSDMData.h"
#include "ObitSource.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSourceEphemerus.h
 *
 * ObitSourceEphemerus template for classes derived from #Obit
 *
 * The ObitSourceEphemerus has position information of moving sources.
 * 
 * \section ObitSourceEphemerusaccess Creators and Destructors
 * An ObitSourceEphemerus will usually be created using ObitSourceEphemerusCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitSourceEphemerus should always be made using the
 * #ObitSourceEphemerusRef function which updates the reference count in the object.
 * Then whenever freeing an ObitSourceEphemerus or changing a pointer, the function
 * #ObitSourceEphemerusUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitSourceEphemerus Class structure. */
typedef struct {
#include "ObitSourceEphemerusDef.h"   /* this class definition */
} ObitSourceEphemerus;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitSourceEphemerus
 * returns a ObitSourceEphemerus*.
 * in = object to unreference
 */
#define ObitSourceEphemerusUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitSourceEphemerus.
 * returns a ObitSourceEphemerus*.
 * in = object to reference
 */
#define ObitSourceEphemerusRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitSourceEphemerusIsA(in) ObitIsA (in, ObitSourceEphemerusGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitSourceEphemerusClassInit (void);

/** Public: Default Constructor. */
ObitSourceEphemerus* newObitSourceEphemerus (gchar* name);

/** Public: Create/initialize ObitSourceEphemerus structures */
ObitSourceEphemerus* ObitSourceEphemerusCreate (gchar* name);
/** Typedef for definition of class pointer structure */
typedef ObitSourceEphemerus* (*ObitSourceEphemerusCreateFP) (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitSourceEphemerusGetClass (void);

/** Public: Copy (deep) constructor. */
ObitSourceEphemerus* ObitSourceEphemerusCopy  (ObitSourceEphemerus *in, 
					       ObitSourceEphemerus *out, 
					       ObitErr *err);

/** Public: Copy structure. */
void ObitSourceEphemerusClone (ObitSourceEphemerus *in, 
			       ObitSourceEphemerus *out, ObitErr *err);

/** Public: Setup given an ObitASDM. */
void ObitSourceEphemerusSetup (ObitSourceEphemerus *in, ObitSDMData *SDM,
			       odouble updtime, ObitUVDesc *uvDesc,
			       ObitErr *err);

/** Public: Check if given source ID included, if so, get position */
gboolean 
ObitSourceEphemerusCheckSource(ObitSourceEphemerus *in, olong srcID,
			       ofloat time, odouble *RA, odouble *Dec, 
			       odouble *dist, ofloat *uvrot);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitSourceEphemerusClassDef.h"
} ObitSourceEphemerusClassInfo; 

#endif /* OBITSOURCEEPHEMERUS_H  */ 
