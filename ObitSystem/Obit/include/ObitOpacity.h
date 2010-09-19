/* $Id: ObitOpacity.h 2 2008-06-10 15:32:27Z bill.cotton $        */
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
#ifndef OBITOPACITY_H 
#define OBITOPACITY_H 

#include "Obit.h"
#include "ObitUV.h"
#include "ObitErr.h"
#include "ObitWeather.h"
#include "ObitSource.h"
#include "ObitSourceList.h"
#include "ObitAntennaList.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOpacity.h
 *
 * ObitOpacity template for classes derived from #Obit
 *
 * The ObitOpacity class calculates atmospheric opacities for a given
 * antenna, frequency, time and direction based on weather or seasonal 
 * information.
 * Uses medhod of VLA Memo # 143 "Improving the frequency resolution 
 * of the default atmospheric opacity model" by Josh Marvil (NRAO), 04/06/2010.
 * 
 * \section ObitOpacityaccess Creators and Destructors
 * An ObitOpacity will usually be created using ObitOpacityCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitOpacity should always be made using the
 * #ObitOpacityRef function which updates the reference count in the object.
 * Then whenever freeing an ObitOpacity or changing a pointer, the function
 * #ObitOpacityUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitOpacity Class structure. */
typedef struct {
#include "ObitOpacityDef.h"   /* this class definition */
} ObitOpacity;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitOpacity
 * returns a ObitOpacity*.
 * in = object to unreference
 */
#define ObitOpacityUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitOpacity.
 * returns a ObitOpacity*.
 * in = object to reference
 */
#define ObitOpacityRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitOpacityIsA(in) ObitIsA (in, ObitOpacityGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitOpacityClassInit (void);

/** Public: Default Constructor. */
ObitOpacity* newObitOpacity (gchar* name);

/** Public: Create/initialize ObitOpacity structures */
ObitOpacity* ObitOpacityCreate (gchar* name, ObitUV *inData);
/** Typedef for definition of class pointer structure */
typedef ObitOpacity* (*ObitOpacityCreateFP) (gchar* name, 
					     ObitUV *inData);
/** Public: Calculate opacities for a set of frequencies */
void ObitOpacityCalc (ObitOpacity *in, ofloat time, olong nfreq, odouble *freqs,
		      olong Ant, olong subA, olong Sou, ofloat *opac, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitOpacityCalcFP) (ObitOpacity *in, ofloat time, 
				   olong nfreq, odouble *freqs,
				   olong Ant, olong subA, olong Sou, 
				   ofloat *opac, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitOpacityGetClass (void);

/** Public: Copy (deep) constructor. */
ObitOpacity* ObitOpacityCopy  (ObitOpacity *in, ObitOpacity *out, ObitErr *err);

/** Public: Copy structure. */
void ObitOpacityClone (ObitOpacity *in, ObitOpacity *out, ObitErr *err);


/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitOpacityClassDef.h"
} ObitOpacityClassInfo; 

#endif /* OBITOPACITY_H */ 
