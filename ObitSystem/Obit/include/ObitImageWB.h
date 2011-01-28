/* $Id$ */
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
#ifndef OBITIMAGEWB_H 
#define OBITIMAGEWB_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitIO.h"
#include "ObitImage.h"
#include "ObitTableCC.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageWB.h
 * ObitImageWB class definition.
 *
 * This class is derived from the #ObitImage class.
 *
 * This class contains the various arrays needed for wideband
 * imaging in the style of Sault & Wieringa 1994 A&AS 108, 585
 * Fit S = S_0 exp(alpha log(nu/nu_0) + beta log(nu/nu_0)*2)
 * using approximation:
 * S = S0(1+ alpha*log(nu/nu_0) + beta * log(nu/nu_0)^2)
 * An ObitImageWB is the front end to persistent disk resident 
 * structures.
 * Access to members is via the member's functions.
 * Both FITS and AIPS cataloged images are supported.
 *
 * \section ObitImageWBaccess Creators and Destructors
 * An ObitImageWB can be created using newObitImageWB which 
 * allows specifying  a name for the object.  
 * This name is used to label messages.
 * The copy constructor ObitImageWBCopy will make a shallow copy
 * of an extant #ObitImageWB.  
 *
 * A copy of a pointer to an ObitImageWB should always be made using the
 * ObitImageWBRef function which updates the reference count in the object.
 * Then whenever freeing an ObitImageWB or changing a pointer, the function
 * ObitImageWBUnref will decrement the reference count and destroy the 
 * object when the reference count hits 0.
 *
 */

/*--------------Class definitions-------------------------------------*/
/** ObitImageWB Class structure. */
typedef struct {
#include "ObitImageWBDef.h"   /* this class definition */
} ObitImageWB;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitImageWB
 * returns a ObitImageWB*.
 * in = object to unreference
 */
#define ObitImageWBUnref(in) ObitUnref ( in)

/** 
 * Macro to reference (update reference count) an ObitImageWB.
 * returns a ObitImageWB*.
 * in = object to reference
 */
#define ObitImageWBRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitImageWBIsA(in) ObitIsA (in, ObitImageWBGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitImageWBClassInit (void);

/** Public: Constructor. */
ObitImageWB* newObitImageWB (gchar* name);

/** Public: Create ImageWB object from description in an ObitInfoList */
ObitImageWB* ObitImageWBFromInfo (gchar *prefix, ObitInfoList *inList, 
				  ObitErr *err);

/** Public: Create ObitImageWB object from ObitImage */
ObitImageWB* ObitImageWBFromImage (ObitImage* in, olong norder, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitImageWBGetClass (void);

/** Public: Copy (shallow) constructor. */
ObitImageWB* 
ObitImageWBCopy  (ObitImageWB *in, ObitImageWB *out, 
		  ObitErr *err);

/** Public: Delete underlying structures. */
ObitImage* ObitImageWBZap  (ObitImage *in, ObitErr *err);

/** Public: Set order of SW spectral imaging */
void ObitImageWBSetOrder (ObitImageWB *image, olong order, ObitErr *err);

/** Public: Create/populate  work arrays for SW spectral imaging. */
void ObitImageWBMakeWork (ObitImageWB *in, ObitErr *err);

/** Public: deallocate work arrays for SW spectral imaging */
void ObitImageWBClearWork (ObitImageWB *image);

/** Public: Decompose image to spectral components. */
void ObitImageWBDecomp (ObitImageWB *in, ObitErr *err);

/** Public: Return image Beam. */
ObitImage* ObitImageWBGetBeam (ObitImage *in, olong beamNo, olong plane[5], 
			       ObitErr *err);

/** Public: Return image Beam order. */
olong ObitImageWBGetBeamOrder (ObitImage *in);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitImageWBClassDef.h"
} ObitImageWBClassInfo; 

#endif /* OBITIMAGEWB_H */ 
