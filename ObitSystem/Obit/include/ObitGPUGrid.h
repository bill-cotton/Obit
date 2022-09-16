/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2021,2022                                          */
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
#ifndef OBITGPUGRID_H 
#define OBITGPUGRID_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitUV.h"
/*#include "ObitImage.h"*/
#include "ObitData.h"
#include "ObitCUDAGrid.h"
#include "ObitCUDAGridInfoDef.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGPUGrid.h
 *
 * ObitGPUGrid GPU enhanced uv gtridding class
 * Uses functions in ObitCUDAGrid
 *
 * Grids are calculated using a GPU.
 * 
 * \section ObitGPUGridaccess Creators and Destructors
 * An ObitGPUGrid will usually be created using ObitGPUGridCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitGPUGrid should always be made using the
 * #ObitGPUGridRef function which updates the reference count in the object.
 * Then whenever freeing an ObitGPUGrid or changing a pointer, the function
 * #ObitGPUGridUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitGPUGrid Class structure. */
typedef struct {
#include "ObitGPUGridDef.h"   /* this class definition */
} ObitGPUGrid;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitGPUGrid
 * returns a ObitGPUGrid*.
 * in = object to unreference
 */
#define ObitGPUGridUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitGPUGrid.
 * returns a ObitGPUGrid*.
 * in = object to reference
 */
#define ObitGPUGridRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGPUGridIsA(in) ObitIsA (in, ObitGPUGridGetClass())

/**
 * Maximum number of vis per read for GPU based Gridding
 */
#define GPU_NVISPIO 8192

/*---------------------------------- Structures ----------------------------*/
/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitGPUGridClassInit (void);

/** Public: Default Constructor. */
ObitGPUGrid* newObitGPUGrid (gchar* name);

/** Public: Create/initialize CUDA GPUGrid structures */
ObitGPUGrid* ObitGPUGridCreate (gchar* name, olong nfacet, Obit *image, ObitUV *UVin);
/** Typedef for definition of class pointer structure */
typedef ObitGPUGrid* (*ObitGPUGridCreateFP) (gchar* name, 
					     olong nfacet, Obit *image, ObitUV *UVin);

/** Public: ClassInfo pointer */
gconstpointer ObitGPUGridGetClass (void);

/** Public: Copy (deep) constructor. */
ObitGPUGrid* ObitGPUGridCopy  (ObitGPUGrid *in, ObitGPUGrid *out, ObitErr *err);

/** Public: Copy structure. */
void ObitGPUGridClone (ObitGPUGrid *in, ObitGPUGrid *out, ObitErr *err);

/* Public: Initialize GPU */
void ObitGPUGridInitGPU (ObitGPUGrid *in, 
			 olong nfacet, olong ifacet, olong nplane,
			 ObitUV *uvdata, ObitErr *err);
typedef void (*ObitGPUGridInitGPUFP) (ObitGPUGrid *in, 
				      olong nfacet, olong ifacet, olong nplane,
				      ObitUV *uvdata, ObitErr *err);

/** Public: copy griding parameters to GPU(CUDA) structues */
void ObitGPUGridSetGPUStruct (ObitGPUGrid *in, Obit *uvgrid, gboolean *chDone,
			      olong ifacet, olong nfacet, olong nplane, ObitUV *UVin, 
			      Obit *imagee, gboolean doBeam, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitGPUGridSetGPUStructFP) (ObitGPUGrid *in, Obit *uvgrid, gboolean *chDone,
					   olong ifacet, olong nfacet, olong nplane,
					   ObitUV *UVin, Obit *imagee,
					   gboolean doBeam, ObitErr *err);

/* Public: Setup GPU */
ollong ObitGPUGridSetGPU (olong gpuno);
typedef ollong (*ObitGPUGridSetGPUFP) (olong gpuno);

/* Public: copy UV data to GPU */
void ObitGPUGrid2GPU (ObitGPUGrid *in, ObitUV *uvdata, ObitErr *err);
typedef void (*ObitGPUGrid2GPUFP) (ObitGPUGrid *in, ObitUV *uvdata, ObitErr *err);

/* Public: Grid UV data */
void ObitGPUGridGrid (ObitGPUGrid *in,  ObitErr *err);
typedef void (*ObitGPUGridGridFP) (ObitGPUGrid *in, ObitErr *err);

/* Public: Flip grid */
void ObitGPUGridFlip (ObitGPUGrid *in, ObitErr *err);
typedef void (*ObitGPUGridFlipFP) (ObitGPUGrid *in, ObitErr *err);

/* Public: copy gridded data to host */
void ObitGPUGrid2Host (ObitGPUGrid *in, olong ifacet, ObitCArray **grid, ObitErr *err);
typedef void (*ObitGPUGrid2HostFP) 
    (ObitGPUGrid *in, olong ifacet, ObitCArray grid, ObitErr *err);

/* Public: Shutdown */
void ObitGPUGridShutdown (ObitGPUGrid *in, ObitUV *uvdata, ObitErr *err);
typedef void (*ObitGPUGridShutdownFP) (ObitGPUGrid *in, ObitUV *uvdata, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitGPUGridClassDef.h"
} ObitGPUGridClassInfo; 

#endif /* OBITGPUGRID_H */ 
