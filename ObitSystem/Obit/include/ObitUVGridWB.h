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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITUVGRIDWB_H 
#define OBITUVGRIDWB_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitUV.h"
#include "ObitFArray.h"
#include "ObitCArray.h"
#include "ObitFFT.h"
#include "ObitImageWB.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVGridWB.h
 * ObitUVGridWB uv data class definition with beams for SW wideband imaging
 *
 * This class is derived from the #ObitUVGrid class.
 *
 * This class is for creating and manipulating a UVGridWB as a memory resident 
 * intermediate entity between ObitUV and ObitImage classes.
 * The beam corresponding to each image should be made first using the
 * same ObitUVGridWB.
 * 
 * \section ObitUVGridWBaccess Creators and Destructors
 * An ObitUVGridWB can be created using newObitUVGridWB which allows specifying 
 * a name for the object.  This name is used to label messages.
 *
 * A copy of a pointer to an ObitUVGridWB should always be made using the
 * #ObitUVGridWBRef function which updates the reference count in the object.
 * Then whenever freeing an ObitUVGridWB or changing a pointer, the function
 * #ObitUVGridWBUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 *
 */

/*--------------Class definitions-------------------------------------*/
/** ObitUVGridWB Class structure. */
typedef struct {
#include "ObitUVGridWBDef.h"   /* this class definition */
} ObitUVGridWB;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitUVGridWB
 * returns a ObitUVGridWB*.
 * in = object to unreference
 */
#define ObitUVGridWBUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitUVGridWB.
 * returns a ObitUVGridWB*.
 * in = object to reference
 */
#define ObitUVGridWBRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitUVGridWBIsA(in) ObitIsA (in, ObitUVGridWBGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitUVGridWBClassInit (void);

/** Public: Constructor. */
ObitUVGridWB* newObitUVGridWB (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitUVGridWBGetClass (void);

/** Public: initialize/reset ObitUVGridWB structures */
void ObitUVGridWBSetup (ObitUVGrid *in, ObitUV *UVin, 
			Obit *image,
			gboolean doBeam, ObitErr *err);

/** Public: Read uv data accumulating to grid */
void ObitUVGridWBReadUV (ObitUVGrid *in, ObitUV *UVin, ObitErr *err);

/** Public: Parallel read uv data accumulating to grid */
void ObitUVGridWBReadUVPar (olong nPar, ObitUVGrid **in, ObitUV **UVin, 
			    ObitErr *err);

/** Public: FFT grid to image plane with gridding correction. */
void ObitUVGridWBFFT2Im (ObitUVGrid *in, Obit *out, ObitErr *err);

/** Public: Parallel FFT grid to image plane with gridding correction. */
void ObitUVGridWBFFT2ImPar (olong nPar, ObitUVGrid **in, Obit **out, 
			    ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitUVGridWBClassDef.h"
} ObitUVGridWBClassInfo; 

#endif /* OBITUVGRIDWB_H */ 
