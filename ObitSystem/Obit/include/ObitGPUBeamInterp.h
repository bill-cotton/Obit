/* $Id:  $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2026                                               */
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
#ifndef OBITGPUBEAMINTERP_H 
#define OBITGPUBEAMINTERP_H 

//#include "Obit.h"
//#include "ObitErr.h"
//#include "ObitUV.h"
//#include "ObitFArray.h"
#include "ObitGPUSkyGeom.h"
#include "ObitGPUFArray.h"
#include "CUDABeamInterp.h"
#include "CUDAFArray.h"
#define STREAM_COUNT 4  /* Number of GPU Streams */
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGPUBeamInterp.h
 *
 * ObitGPUBeamInterp GPU enhanced FArray interpolation class
 * Uses functions in ObitCUDAFInterpolate
 *
 * 2D FArrays are interpolated using a GPU.
 * Version for interpolating antenna beam images
 * Coordinates with CUDABeamInterp.
 * 
 * \section ObitGPUBeamInterpaccess Creators and Destructors
 * An ObitGPUBeamInterp will usually be created using ObitGPUBeamInterpCreate which allows 
 * specifying a name for the object as well as other information.
 * The usable member data is in the  CUDAFInterpolate, beamInterp member.
 *
 * A copy of a pointer to an ObitGPUBeamInterp should always be made using the
 * #ObitGPUBeamInterpRef function which updates the reference count in the object.
 * Then whenever freeing an ObitGPUBeamInterp or changing a pointer, the function
 * #ObitGPUBeamInterpUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
#if HAVE_GPU==1  /* Have a GPU? */
#ifndef CUDABEAMINTERPDEF_H // Prevent multiple definitions
#define CUDABEAMINTERPDEF_H
typedef struct {
#include "CUDABeamInterpDef.h"   /* CUDA stuff */
} CUDABeamInterp;
#endif /* CUDABEAMINTERPDEF_H */
#endif /* HAVE_GPU */

/** ObitGPUBeamInterp Class structure. */
typedef struct {
#include "ObitGPUBeamInterpDef.h"   /* this class definition */
} ObitGPUBeamInterp;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitGPUBeamInterp
 * returns a ObitGPUBeamInterp*.
 * in = object to unreference
 */
#define ObitGPUBeamInterpUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitGPUBeamInterp.
 * returns a ObitGPUBeamInterp*.
 * in = object to reference
 */
#define ObitGPUBeamInterpRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGPUBeamInterpIsA(in) ObitIsA (in, ObitGPUBeamInterpGetClass())

/*---------------------------------- Structures ----------------------------*/
/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitGPUBeamInterpClassInit (void);

/** Public: Default Constructor. */
ObitGPUBeamInterp* newObitGPUBeamInterp (gchar* name);

/** Public: Create/initialize ObitGPUBeamInterp structures */
ObitGPUBeamInterp* ObitGPUBeamInterpCreate (gchar* name, 
					    ObitFArray *inArray,  ObitImageDesc *Desc, 
					    int hwidth, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitGPUBeamInterp* 
(*ObitGPUBeamInterpCreateFP) (gchar* name,
			      ObitFArray *inArray,  ObitImageDesc *Desc, 
			      int hwidth, ObitErr *err);

#if HAVE_GPU==1  /* Have a GPU? */
/* Public: Create CUDABeamInterp and copy parts from host */
CUDABeamInterp* ObitGPUCUDABeamInterpCreate (ObitFArray *array, ObitImageDesc *Desc,
					     int hwidth);
/** Typedef for definition of class pointer structure */
typedef CUDABeamInterp*  (*ObitGPUCUDABeamInterpCreateFP)
  (ObitFArray *array, ObitImageDesc *Desc, int hwidth);

/** Free a CUDABeamInterp */
void CUDABeamInterpFree (CUDABeamInterp* in);

#else /* not GPU enabled */
/* Public: Create CUDABeamInterp and copy parts from host */
void* ObitGPUCUDABeamInterpCreate (ObitFArray *array, ObitImageDesc *Desc,
				   int hwidth);
/** Typedef for definition of class pointer structure */
typedef void*  (*ObitGPUCUDABeamInterpCreateFP)
  (ObitFArray *array, ObitImageDesc *Desc, int hwidth);

/** Free a CUDABeamInterp */
void CUDABeamInterpFree (void* in);

#endif /* HAVE_GPU */
/** Public: Host to Device Image Descriptor conversion */
CUDAImageDesc* CUDASkyGeomImageH2D (ObitImageDesc *in);

/** Public: Host to Device UV Descriptor conversion */
CUDAUVDesc* CUDASkyGeomUVH2D (ObitUVDesc *in);

/** Public: ClassInfo pointer */
gconstpointer ObitGPUBeamInterpGetClass (void);

/** Public: Copy (deep) constructor. */
ObitGPUBeamInterp* ObitGPUBeamInterpCopy  (ObitGPUBeamInterp *in, 
					       ObitGPUBeamInterp *out, ObitErr *err);

/** Public: Copy structure. */
void ObitGPUBeamInterpClone (ObitGPUBeamInterp *in, 
			       ObitGPUBeamInterp *out, ObitErr *err);


/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitGPUBeamInterpClassDef.h"
} ObitGPUBeamInterpClassInfo; 

#endif /* OBITGPUBEAMINTERP_H  */ 
