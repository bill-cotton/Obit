/* $Id: ObitWeather.h 2 2008-06-10 15:32:27Z bill.cotton $        */
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
#ifndef OBITWEATHER_H 
#define OBITWEATHER_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitUV.h"
#include "ObitWeather.h"
#include "ObitTableWX.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitWeather.h
 *
 * ObitWeather Weather table interpolator class
 *
 * The ObitWeather class provides an interface to a weather table allowing
 * interpolation for a given time and antenna.
 * 
 * \section ObitWeatheraccess Creators and Destructors
 * An ObitWeather will usually be created using ObitWeatherCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitWeather should always be made using the
 * #ObitWeatherRef function which updates the reference count in the object.
 * Then whenever freeing an ObitWeather or changing a pointer, the function
 * #ObitWeatherUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitWeather Class structure. */
typedef struct {
#include "ObitWeatherDef.h"   /* this class definition */
} ObitWeather;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitWeather
 * returns a ObitWeather*.
 * in = object to unreference
 */
#define ObitWeatherUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitWeather.
 * returns a ObitWeather*.
 * in = object to reference
 */
#define ObitWeatherRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitWeatherIsA(in) ObitIsA (in, ObitWeatherGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitWeatherClassInit (void);

/** Public: Default Constructor. */
ObitWeather* newObitWeather (gchar* name);

/** Public: Create/initialize ObitWeather structures */
ObitWeather* ObitWeatherCreate (gchar* name, ObitTableWX *WXTable, 
				ObitUV *UVData, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitWeather* (*ObitWeatherCreateFP) (gchar* name, ObitTableWX *WXTable, 
					     ObitUV *UVData, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitWeatherGetClass (void);

/** Public: Copy (deep) constructor. */
ObitWeather* ObitWeatherCopy  (ObitWeather *in, ObitWeather *out, ObitErr *err);

/** Public: Copy structure. */
void ObitWeatherClone (ObitWeather *in, ObitWeather *out, ObitErr *err);

/** Public: Interpolate Weather values. */
void ObitWeatherReport (ObitWeather *in, ofloat time, olong ant, olong suba,
			ofloat *temp, ofloat *DP, ofloat *press, 
			ofloat *windDir, ofloat *windVel, 
			ofloat *wvrH2O, ofloat *ions, 
			ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitWeatherClassDef.h"
} ObitWeatherClassInfo; 

#endif /* OBITWEATHER_H */ 
