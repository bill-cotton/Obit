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
#ifndef OBITDOPPLER_H 
#define OBITDOPPLER_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitFile.h"
#include "ObitThread.h"
#include "ObitUV.h"
#include "ObitFFT.h"
#include "ObitInfoList.h"
#include "ObitAntennaList.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDoppler.h
 *
 * ObitDoppler UV data Doppler correcting class for  #Obit
 *
 * \section ObitDoppleraccess Creators and Destructors
 * An ObitDoppler will usually be created using ObitDopplerCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitDoppler should always be made using the
 * #ObitDopplerRef function which updates the reference count in the object.
 * Then whenever freeing an ObitDoppler or changing a pointer, the function
 * #ObitDopplerUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitDoppler Class structure. */
typedef struct {
#include "ObitDopplerDef.h"   /* this class definition */
} ObitDoppler;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitDoppler
 * returns a ObitDoppler*.
 * in = object to unreference
 */
#define ObitDopplerUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitDoppler.
 * returns a ObitDoppler*.
 * in = object to reference
 */
#define ObitDopplerRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitDopplerIsA(in) ObitIsA (in, ObitDopplerGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitDopplerClassInit (void);

/** Public: Default Constructor. */
ObitDoppler* newObitDoppler (gchar* name);

/** Public: Create/initialize ObitDoppler structures */
ObitDoppler* ObitDopplerCreate (gchar* name, ObitUV *uvdata, 
				ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitDoppler* (*ObitDopplerCreateFP) (gchar* name,  
					     ObitUV *uvdata, 
					     ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitDopplerGetClass (void);

/** Public: Copy (deep) constructor. */
ObitDoppler* ObitDopplerCopy  (ObitDoppler *in, ObitDoppler *out, ObitErr *err);

/** Public: Copy structure. */
void ObitDopplerClone (ObitDoppler *in, ObitDoppler *out, ObitErr *err);

/** Public: Calculate observed frequency from a given line, source, time...  */
odouble ObitDopplerFreqLSR (odouble rest, ofloat vlsr, 
			    odouble ra, odouble dec, 
			    olong year, olong doy, odouble ut, 
			    odouble x, odouble y, odouble z);  

/** Public: Make Doppler corrections to UV data set */
ObitUV* ObitDopplerCVel (ObitUV *inUV, gboolean scratch, ObitUV *outUV, 
			 ObitErr *err);
/** Public: Convert JD to year, doy, time */
void ObitDopplerJD2Date (odouble JD, olong *year, olong *doy, ofloat *ut);
/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitDopplerClassDef.h"
} ObitDopplerClassInfo; 

#endif /* OBITFDOPPLER_H */ 
