/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2010                                          */
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
#ifndef OBITUVGRID_H 
#define OBITUVGRID_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitUV.h"
#include "ObitFArray.h"
#include "ObitCArray.h"
#include "ObitFFT.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVGrid.h
 * ObitUVGrid uv data class definition.
 *
 * This class is derived from the #Obit class.
 *
 * This class is for creating and manipulating a UVGrid as a memory resident 
 * intermediate entity between ObitUV and ObitImage classes.
 * The beam corresponding to each image should be made first using the
 * same ObitUVGrid.
 * 
 * \section ObitUVGridaccess Creators and Destructors
 * An ObitUVGrid can be created using newObitUVGrid which allows specifying 
 * a name for the object.  This name is used to label messages.
 *
 * A copy of a pointer to an ObitUVGrid should always be made using the
 * #ObitUVGridRef function which updates the reference count in the object.
 * Then whenever freeing an ObitUVGrid or changing a pointer, the function
 * #ObitUVGridUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 *
 */

/*--------------Class definitions-------------------------------------*/
/** ObitUVGrid Class structure. */
typedef struct {
#include "ObitUVGridDef.h"   /* this class definition */
} ObitUVGrid;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitUVGrid
 * returns a ObitUVGrid*.
 * in = object to unreference
 */
#define ObitUVGridUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitUVGrid.
 * returns a ObitUVGrid*.
 * in = object to reference
 */
#define ObitUVGridRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitUVGridIsA(in) ObitIsA (in, ObitUVGridGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitUVGridClassInit (void);

/** Public: Constructor. */
ObitUVGrid* newObitUVGrid (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitUVGridGetClass (void);

/** Public: initialize/reset ObitUVGrid structures */
void ObitUVGridSetup (ObitUVGrid *in, ObitUV *UVin, 
		       Obit *image,
		       gboolean doBeam, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitUVGridSetupFP) (ObitUVGrid *in, ObitUV *UVin, 
		       Obit *image,
		       gboolean doBeam, ObitErr *err);

/** Public: Read uv data accumulating to grid */
void ObitUVGridReadUV (ObitUVGrid *in, ObitUV *UVin, ObitErr *err);
typedef ObitIOCode (*ObitUVGridReadUVFP) (ObitUVGrid *in, ObitUV *UVin, 
					 ObitErr *err);

/** Public: Parallel read uv data accumulating to grid */
void ObitUVGridReadUVPar (olong nPar, ObitUVGrid **in, ObitUV **UVin, ObitErr *err);
typedef ObitIOCode (*ObitUVGridReadUVParFP) (olong nPar, ObitUVGrid **in, ObitUV **UVin, 
					     ObitErr *err);

/** Public: FFT grid to image plane with gridding correction. */
void ObitUVGridFFT2Im (ObitUVGrid *in, Obit *out, ObitErr *err);
typedef void (*ObitUVGridFFT2ImFP) (ObitUVGrid *in, Obit *out, ObitErr *err);

/** Public: Parallel FFT grid to image plane with gridding correction. */
void ObitUVGridFFT2ImPar (olong nPar, ObitUVGrid **in, Obit **out, ObitErr *err);
typedef void (*ObitUVGridFFT2ImParFP) (olong nPar, ObitUVGrid **in, Obit **out, 
				       ObitErr *err);

/* Private functions for derived classes */
/** Prepare uv buffer for gridding */
void PrepBuffer (ObitUVGrid* in, ObitUV *uvdata, olong loVis, olong hiVis);
typedef void (*PrepBufferFP) (ObitUVGrid* in, ObitUV *uvdata, olong loVis, olong hiVis);

/** Grid uv buffer */
void GridBuffer (ObitUVGrid* in, ObitUV *uvdata, olong loVis, olong hiVis,
		 ObitCArray *accGrid);
typedef void (*GridBufferFP) (ObitUVGrid* in, ObitUV *uvdata, olong loVis, olong hiVis,
	      ObitCArray *accGrid);

/** Gridding correction */
void GridCorrFn (ObitUVGrid* in, long n, olong icent, 
		 ofloat *data, ofloat *ramp, ObitFArray *out);
typedef void (*GridCorrFnFP) (ObitUVGrid* in, long n, olong icent, 
			      ofloat *data, ofloat *ramp, ObitFArray *out);
/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitUVGridClassDef.h"
} ObitUVGridClassInfo; 

#endif /* OBITUVGRID_H */ 
