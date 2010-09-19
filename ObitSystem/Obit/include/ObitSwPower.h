/* $Id: ObitSwPower.h 2 2008-06-10 15:32:27Z bill.cotton $        */
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
#ifndef OBITSWPOWER_H 
#define OBITSWPOWER_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitUV.h"
#include "ObitSwPower.h"
#include "ObitTableSY.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSwPower.h
 *
 * ObitSwPower EVLA Switched Power interpolator class
 *
 * The ObitSwPower class provides an interface to a EVLA Switched power (SY)
 * table allowing interpolation for a given time and antenna.
 * 
 * \section ObitSwPoweraccess Creators and Destructors
 * An ObitSwPower will usually be created using ObitSwPowerCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitSwPower should always be made using the
 * #ObitSwPowerRef function which updates the reference count in the object.
 * Then whenever freeing an ObitSwPower or changing a pointer, the function
 * #ObitSwPowerUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitSwPower Class structure. */
typedef struct {
#include "ObitSwPowerDef.h"   /* this class definition */
} ObitSwPower;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitSwPower
 * returns a ObitSwPower*.
 * in = object to unreference
 */
#define ObitSwPowerUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitSwPower.
 * returns a ObitSwPower*.
 * in = object to reference
 */
#define ObitSwPowerRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitSwPowerIsA(in) ObitIsA (in, ObitSwPowerGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitSwPowerClassInit (void);

/** Public: Default Constructor. */
ObitSwPower* newObitSwPower (gchar* name);

/** Public: Create/initialize ObitSwPower structures */
ObitSwPower* ObitSwPowerCreate (gchar* name, ObitTableSY *SYTable, 
			  ObitUV *UVData, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitSwPower* (*ObitSwPowerCreateFP) (gchar* name, ObitTableSY *SYTable, 
				       ObitUV *UVData, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitSwPowerGetClass (void);

/** Public: Copy (deep) constructor. */
ObitSwPower* ObitSwPowerCopy  (ObitSwPower *in, ObitSwPower *out, ObitErr *err);

/** Public: Copy structure. */
void ObitSwPowerClone (ObitSwPower *in, ObitSwPower *out, ObitErr *err);

/** Public: Interpolate SwPower values. */
void ObitSwPowerReport (ObitSwPower *in, ofloat time, olong ant, olong suba,
			ofloat *PwrDif1, ofloat *PwrSum1, ofloat *Gain1,
			ofloat *PwrDif2, ofloat *PwrSum2, ofloat *Gain2,
			ObitErr *err);

/** Public: Smooth SY table. */
void ObitSwPowerSYSmo (ObitTableSY *SYTab, olong isuba, ObitErr* err) ;

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitSwPowerClassDef.h"
} ObitSwPowerClassInfo; 

#endif /* OBITWPOWER_H */ 
