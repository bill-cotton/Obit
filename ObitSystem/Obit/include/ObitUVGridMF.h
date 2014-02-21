/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010-2014                                          */
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
#ifndef OBITUVGRIDMF_H 
#define OBITUVGRIDMF_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitUV.h"
#include "ObitFArray.h"
#include "ObitCArray.h"
#include "ObitFFT.h"
#include "ObitImageMF.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVGridMF.h
 * ObitUVGridMF uv data class definition with beams for wideband spectral imaging
 *
 * This class is derived from the #ObitUVGrid class.
 *
 * This class is for creating and manipulating a UVGridMF as a memory resident 
 * intermediate entity between ObitUV and ObitImage classes.
 * The beam corresponding to each image should be made first using the
 * same ObitUVGridMF.
 * 
 * \section ObitUVGridMFaccess Creators and Destructors
 * An ObitUVGridMF can be created using newObitUVGridMF which allows specifying 
 * a name for the object.  This name is used to label messages.
 *
 * A copy of a pointer to an ObitUVGridMF should always be made using the
 * #ObitUVGridMFRef function which updates the reference count in the object.
 * Then whenever freeing an ObitUVGridMF or changing a pointer, the function
 * #ObitUVGridMFUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 *
 */

/*--------------Class definitions-------------------------------------*/
/** ObitUVGridMF Class structure. */
typedef struct {
#include "ObitUVGridMFDef.h"   /* this class definition */
} ObitUVGridMF;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitUVGridMF
 * returns a ObitUVGridMF*.
 * in = object to unreference
 */
#define ObitUVGridMFUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitUVGridMF.
 * returns a ObitUVGridMF*.
 * in = object to reference
 */
#define ObitUVGridMFRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitUVGridMFIsA(in) ObitIsA (in, ObitUVGridMFGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitUVGridMFClassInit (void);

/** Public: Constructor. */
ObitUVGridMF* newObitUVGridMF (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitUVGridMFGetClass (void);

/** Public: initialize/reset ObitUVGridMF structures */
void ObitUVGridMFSetup (ObitUVGrid *in, ObitUV *UVin, 
			Obit *image,
			gboolean doBeam, ObitErr *err);

/** Public: Read uv data accumulating to grid */
void ObitUVGridMFReadUV (ObitUVGrid *in, ObitUV *UVin, ObitErr *err);

/** Public: Parallel read uv data accumulating to grid */
void ObitUVGridMFReadUVPar (olong nPar, ObitUVGrid **in, ObitUV *UVin, 
			    ObitErr *err);

/** Public: FFT grid to image plane with gridding correction. */
void ObitUVGridMFFFT2Im (ObitUVGrid *in, Obit *out, ObitErr *err);

/** Public: Parallel FFT grid to image plane with gridding correction. */
void ObitUVGridMFFFT2ImPar (olong nPar, ObitUVGrid **in, Obit **out, 
			    ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitUVGridMFClassDef.h"
} ObitUVGridMFClassInfo; 

#endif /* OBITUVGRIDMF_H */ 
