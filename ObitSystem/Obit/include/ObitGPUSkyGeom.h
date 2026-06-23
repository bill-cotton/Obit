/* $Id: $   */

/*HIDE c (esp.glib) structures from cuda */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2026                                               */
/*;  Associated Universities, Inc. Washington DC, USA.                */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the ObitGPUSkyGeom (GPU version of ObitImageDesc+ObitSkyGeom) */
/**
 * \file ObitGPUSkyGeom.h
 * GPU version of ObitImageDesc+ObitSkyGeom utilities
 */
#ifndef  OBITGPUSKYGEOM_H 
#define  OBITGPUSKYGEOM_H
#include "ObitImageDesc.h"
#include "ObitUVDesc.h"
#include "ObitFArray.h"

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitCoordType
 * enum for coordinate system types (Equatorial, Galactic, Ecliptic)
 */
enum cudaCoordType {
  /** Equatorial (RA, Dec) */
  CUDA_Equatorial, 
  /** Galactic (GLONG, GLAT) */
  CUDA_Galactic,  
  /** Ecliptic (ELONG, ELAT) */
  CUDA_Ecliptic
}; /* end enum cudaCoordType */
/** typedef for enum for coordinate system types. */
typedef enum cudaCoordType CUDACoordType;

/*--------------Class definitions-------------------------------------*/
/** ObitGPUSkyGeom Class structures: CUDAImageDesc, CUDAUVDesc */
#include "ObitGPUSkyGeomDef.h"   /* this class definition */

/** Public: Create/initialize CUDAImageDesc structures */
CUDAImageDesc* ObitGPUImageDescCreate ();

/** Public: ImageDesc Destructor */
void ObitGPUImageDescZap (CUDAImageDesc *in);

/** Public: Create/initialize CUDAUVDesc structures */
CUDAUVDesc* ObitGPUUVDescCreate ();

/** Public: UVDesc Destructor */
void ObitGPUVDescZap (CUDAUVDesc *in);

/** Public: Host to Device Image Descriptor conversion */
CUDAImageDesc* ObitGPUSkyGeomImageH2D (ObitImageDesc *in);

/** Public: Host to Device UV Descriptor conversion */
CUDAUVDesc* ObitGPUSkyGeomUVH2D (ObitUVDesc *in);

#endif  /* OBITGPUSKYGEOM_H */
