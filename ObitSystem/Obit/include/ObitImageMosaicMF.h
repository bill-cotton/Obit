/* $Id$ */
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
#ifndef OBITIMAGEMOSAICMF_H 
#define OBITIMAGEMOSAICMF_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitIO.h"
#include "ObitImage.h"
#include "ObitTableCC.h"
#include "ObitImageMosaic.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageMosaicMF.h
 * ObitImageMosaicMF class definition.
 *
 * This class is derived from the #ObitImageMosaic class.
 *
 * An ObitImageMosaicMF is an array of associated images, generally 
 * intended to cover a region of the sky.
 * This class contains the various arrays needed for wideband
 * imaging in the spectral imaging style
 * Fit S = S_0 exp(alpha log(nu/nu_0) + beta log(nu/nu_0)*2)
 * This class contains an array of astronomical images and allows access.
 * An ObitImageMosaicMF is the front end to persistent disk resident 
 * structures.
 * Access to members is via the member's functions.
 * Both FITS and AIPS cataloged images are supported.
 *
 * \section ObitImageMosaicMFaccess Creators and Destructors
 * An ObitImageMosaicMF can be created using newObitImageMosaicMF which 
 * allows specifying  a name for the object.  
 * This name is used to label messages.
 * The copy constructor ObitImageMosaicMFCopy will make a shallow copy
 * of an extant #ObitImageMosaicMF.  
 *
 * A copy of a pointer to an ObitImageMosaicMF should always be made using the
 * ObitImageMosaicMFRef function which updates the reference count in the object.
 * Then whenever freeing an ObitImageMosaicMF or changing a pointer, the function
 * ObitImageMosaicMFUnref will decrement the reference count and destroy the 
 * object when the reference count hits 0.
 *
 */

/*--------------Class definitions-------------------------------------*/
/** ObitImageMosaicMF Class structure. */
typedef struct {
#include "ObitImageMosaicMFDef.h"   /* this class definition */
} ObitImageMosaicMF;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitImageMosaicMF
 * returns a ObitImageMosaicMF*.
 * in = object to unreference
 */
#define ObitImageMosaicMFUnref(in) ObitUnref ( in)

/** 
 * Macro to reference (update reference count) an ObitImageMosaicMF.
 * returns a ObitImageMosaicMF*.
 * in = object to reference
 */
#define ObitImageMosaicMFRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitImageMosaicMFIsA(in) ObitIsA (in, ObitImageMosaicMFGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitImageMosaicMFClassInit (void);

/** Public: Constructor. */
ObitImageMosaicMF* newObitImageMosaicMF (gchar* name, olong number);

/** Public: Create ImageMosaicMF object from description in an ObitInfoList */
ObitImageMosaicMF* ObitImageMosaicMFFromInfo (gchar *prefix, ObitInfoList *inList, 
					      ObitErr *err);
typedef ObitImageMosaicMF* (*ObitImageMosaicMFFromInfoFP) (gchar *prefix, 
							   ObitInfoList *inList, 
							   ObitErr *err);
/** Public: ClassInfo pointer */
gconstpointer ObitImageMosaicMFGetClass (void);

/** Public: Copy (shallow) constructor. */
ObitImageMosaic* 
ObitImageMosaicMFCopy  (ObitImageMosaic *in, ObitImageMosaic *out, 
			ObitErr *err);

/** Public: Zap specified image. */
void 
ObitImageMosaicMFZapImage  (ObitImageMosaic *in, olong number,
			    ObitErr *err);

/** Public: Set underlying files */
void ObitImageMosaicMFSetFiles  (ObitImageMosaic *in, gboolean doBeam, ObitErr *err);

/** Public: Create Mosaic from uv data */
ObitImageMosaicMF *ObitImageMosaicMFCreate (gchar *name, olong order, ofloat maxFBW,
					    ofloat alpha, odouble alphaRefF,
					    ObitUV *uvData, ObitErr *err);

/** Public: Define parameters of images */
void ObitImageMosaicMFDefine (ObitImageMosaic *in, ObitUV *uvData, gboolean doBeam,
			      ObitErr *err);

/** Public: Flatten tiles onto full field image */
void ObitImageMosaicMFFlatten (ObitImageMosaic *in, ObitErr *err);

/** Public: Reimaging needed to center strong source on pixel? */
gboolean ObitImageMosaicMFReimage (ObitImageMosaicMF *mosaic, ObitErr* err);

/** Public: Add field to mosaic */
void ObitImageMosaicMFAddField (ObitImageMosaic *in, ObitUV *uvData, 
				olong nx, olong ny, olong nplane, 
				ofloat RAShift, ofloat DecShift, 
				gboolean isAuto, ObitErr *err);

/** Public:  Generate a mosaic for peeling */
ObitImageMosaicMF* ObitImageMosaicMFMaxField (ObitImageMosaic *mosaic, 
					      ofloat MinFlux, olong *ignore, olong *field,
					      ObitErr* err); 

/** Public: Extract information about underlying structures to ObitInfoList */
void ObitImageMosaicMFGetInfo (ObitImageMosaic *in, gchar *prefix, ObitInfoList *outList, 
			       ObitErr *err);

/** Public: Concatenate Image CC tables onto the FullField Image */
void ObitImageMosaicMFCopyCC (ObitImageMosaic *in, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitImageMosaicMFClassDef.h"
} ObitImageMosaicMFClassInfo; 

#endif /* OBITIMAGEMOSAICMF_H */ 
