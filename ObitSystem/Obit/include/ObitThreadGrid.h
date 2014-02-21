/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014                                               */
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
#ifndef OBITTHREADGRID_H 
#define OBITTHREADGRID_H 
#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitUV.h"
#include "ObitFArray.h"
#include "ObitCArray.h"
#include "ObitUVGrid.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitThreadGrid.h
 * ObitThreadGrid uv data threaded gridding class definition.
 *
 * This class is derived from the #Obit class.
 *
 * This class is for gridding uv data possibly using threads and/or AVX.
 * 
 * \section ObitThreadGridaccess Creators and Destructors
 * An ObitThreadGrid can be created using newObitThreadGrid which allows specifying 
 * a name for the object.  This name is used to label messages.
 *
 * A copy of a pointer to an ObitThreadGrid should always be made using the
 * #ObitThreadGridRef function which updates the reference count in the object.
 * Then whenever freeing an ObitThreadGrid or changing a pointer, the function
 * #ObitThreadGridUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 *
 */

/*--------------Class definitions-------------------------------------*/
/** ObitThreadGrid Class structure. */
/* Gridding info for */
typedef struct  {
/** Number of facets */
olong nfacet;
/** Number of threads */
olong nThreads;
/** Number of channels */
olong nchan;
/** Number of random parameters in visibilities */
olong nrparm;
/** Number of visibilities */
olong nvis;
/** length in floats of visibilities */
olong lenvis;
/** Convolution kernal width */
olong convWidth;
/** Number of convolution entries per cell */
olong convNperCell;
/** guardband in u, v in lambda per facet */
ofloat *guardu, *guardv;
/** max, min BL per facet */
gfloat *maxBL, *minBL;
/** Size of facets (cells) per facet */
olong *nx, *ny;
/** Size of each grid (nx x ny/2+1 x 2 floats) per facet */
olong *sizeGrid;
/** Scaling for u,v to cells per facet */
ofloat *uscale, *vscale;
/** Circular Beam taper kLambda per facet */
gfloat *BeamTaperUV;
/** shift parameters per facet [3] */
ofloat *shift;
/** UVW, Rotation matrix per facet [3][3] */
ofloat *rotUV;
/** Frequency scaling array (nchan) */
ofloat *freqArr;
/** gridding kernal ordered by fractional cell then cell*/
ofloat *convfn;
/** If facet a beam, per facet */
gboolean *isBeam;
/** Thread to use */
ObitThread *thread;
/** Array of thread arguments */
gpointer *thArgs;
} ObitThreadGridInfo; 

/*--------------Class definitions-------------------------------------*/
/** ObitThreadGrid Class structure. */
typedef struct {
#include "ObitThreadGridDef.h"   /* this class definition */
} ObitThreadGrid;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitThreadGrid
 * returns a ObitThreadGrid*.
 * in = object to unreference
 */
#define ObitThreadGridUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitThreadGrid.
 * returns a ObitThreadGrid*.
 * in = object to reference
 */
#define ObitThreadGridRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitThreadGridIsA(in) ObitIsA (in, ObitThreadGridGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitThreadGridClassInit (void);

/** Public: Constructor. */
ObitThreadGrid* newObitThreadGrid (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitThreadGridGetClass (void);

/** Public: initialize/reset ObitThreadGrid structures */
void ObitThreadGridSetup (ObitThreadGrid *in, ObitUV *UVin,
			  olong nPar, ObitUVGrid **UVGrids, 
			  olong nThreads, ObitErr *err);

/** Typedef for definition of class pointer structure */
typedef void (*ObitThreadGridSetupFP) (ObitThreadGrid *in, ObitUV *UVin,
				       olong nPar, ObitUVGrid **UVGrids, 
				       olong nThreads, ObitErr *err);

/** Public: Initialize basic object */
void ObitThreadGridInit (gpointer out);
/** Public: grid a visibility buffer possible with threads */
void ObitThreadGridGrid (ObitThreadGrid *in);
typedef  void (*ObitThreadGridGridFP)  (ObitThreadGrid *in);
/** Public: flip/add conjugate columns of a set of grids */
void ObitThreadGridFlip (ObitThreadGrid *in);
typedef  void (*ObitThreadGridFlipFP) (ObitThreadGrid *in);
/** Public: swap/Merge versions of a grid from separate threads */
void ObitThreadGridMerge (ObitThreadGrid *in);
typedef  void (*ObitThreadGridMergeFP) (ObitThreadGrid *in);
/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitThreadGridClassDef.h"
} ObitThreadGridClassInfo; 

#endif /* OBITTHREADGRID_H */ 
