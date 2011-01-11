/* $Id: ObitImageMF.h 128 2009-09-23 14:48:29Z bill.cotton $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010,2011                                          */
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
#ifndef OBITIMAGEMF_H 
#define OBITIMAGEMF_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitIO.h"
#include "ObitImage.h"
#include "ObitTableCC.h"
#include "ObitUV.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageMF.h
 * ObitImageMF class definition.
 *
 * This class is derived from the #ObitImage class.
 *
 * This class contains the various arrays needed for wideband
 * imaging in the spectral imaging style.
 * Fit S = S_0 exp(alpha log(nu/nu_0) + beta log(nu/nu_0)*2)
 * An ObitImageMF is the front end to persistent disk resident 
 * structures.
 * Access to members is via the member's functions.
 * Both FITS and AIPS cataloged images are supported.
 *
 * \section ObitImageMFaccess Creators and Destructors
 * An ObitImageMF can be created using newObitImageMF which 
 * allows specifying  a name for the object.  
 * This name is used to label messages.
 * The copy constructor ObitImageMFCopy will make a shallow copy
 * of an extant #ObitImageMF.  
 *
 * A copy of a pointer to an ObitImageMF should always be made using the
 * ObitImageMFRef function which updates the reference count in the object.
 * Then whenever freeing an ObitImageMF or changing a pointer, the function
 * ObitImageMFUnref will decrement the reference count and destroy the 
 * object when the reference count hits 0.
 *
 */

/*--------------Class definitions-------------------------------------*/
/** ObitImageMF Class structure. */
typedef struct {
#include "ObitImageMFDef.h"   /* this class definition */
} ObitImageMF;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitImageMF
 * returns a ObitImageMF*.
 * in = object to unreference
 */
#define ObitImageMFUnref(in) ObitUnref ( in)

/** 
 * Macro to reference (update reference count) an ObitImageMF.
 * returns a ObitImageMF*.
 * in = object to reference
 */
#define ObitImageMFRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitImageMFIsA(in) ObitIsA (in, ObitImageMFGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitImageMFClassInit (void);

/** Public: Constructor. */
ObitImageMF* newObitImageMF (gchar* name);

/** Public: Create ImageMF object from description in an ObitInfoList */
ObitImageMF* ObitImageMFFromInfo (gchar *prefix, ObitInfoList *inList, 
				  ObitErr *err);

/** Public: Create ObitImageMF object from ObitImage */
ObitImageMF* ObitImageMFFromImage (ObitImage* in, ObitUV *inData,
				   olong norder, ofloat maxFBW, 
				   ofloat alpha, odouble alphaRefF,
				   ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitImageMFGetClass (void);

/** Public: Copy (shallow) constructor. */
ObitImageMF* 
ObitImageMFCopy  (ObitImageMF *in, ObitImageMF *out, 
		  ObitErr *err);

/** Public: Delete underlying structures. */
ObitImage* ObitImageMFZap  (ObitImage *in, ObitErr *err);

/** Public: Set Coarse spectral channels */
void ObitImageMFSetSpec (ObitImageMF *image, ObitUV *inData, ofloat maxFBW,
			 ofloat alpha, odouble alphaRefF, ObitErr *err);

/** Public: Get Coarse spectral channel info */
void ObitImageMFGetSpec (ObitImageMF *image, ObitErr *err);

/** Public: Set order of SW spectral imaging */
void ObitImageMFSetOrder (ObitImageMF *image, olong order, ObitErr *err);

/** Public: Make image or beam from coarse channels. */
void ObitImageMFCombine (ObitImageMF *in, gboolean addExt, ObitErr *err);

/** Public: Blanks combined image and spectral terms  */
void ObitImageMFBlank (ObitImageMF *in, ObitErr *err);

/** Public: Fit coarse spectrum planes and write spectral  */
void ObitImageMFFitSpec (ObitImageMF *in, ofloat antSize, ObitErr *err);

/** Public: Return image Beam. */
ObitImage* ObitImageMFGetBeam (ObitImage *in, olong beamNo, olong plane[5], 
			       ObitErr *err);

/** Public: Return image Beam order. */
olong ObitImageMFGetBeamOrder (ObitImage *in);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitImageMFClassDef.h"
} ObitImageMFClassInfo; 

#endif /* OBITIMAGEMF_H */ 
