/* $Id: ObitImageMosaicWB.h 128 2009-09-23 14:48:29Z bill.cotton $ */
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
#ifndef OBITIMAGEMOSAICWB_H 
#define OBITIMAGEMOSAICWB_H 

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
 * \file ObitImageMosaicWB.h
 * ObitImageMosaicWB class definition.
 *
 * This class is derived from the #ObitImageMosaic class.
 *
 * An ObitImageMosaicWB is an array of associated images, generally 
 * intended to cover a region of the sky.
 * This class contains the various arrays needed for wideband
 * imaging in the style of Sault & Wieringa 1994 A&AS 108, 585
 * Fit S = S_0 exp(alpha log(nu/nu_0) + beta log(nu/nu_0)*2)
 * using approximation:
 * S = S0(1+ alpha*log(nu/nu_0) + beta * log(nu/nu_0)^2)
 * This class contains an array of astronomical images and allows access.
 * An ObitImageMosaicWB is the front end to persistent disk resident 
 * structures.
 * Access to members is via the member's functions.
 * Both FITS and AIPS cataloged images are supported.
 *
 * \section ObitImageMosaicWBaccess Creators and Destructors
 * An ObitImageMosaicWB can be created using newObitImageMosaicWB which 
 * allows specifying  a name for the object.  
 * This name is used to label messages.
 * The copy constructor ObitImageMosaicWBCopy will make a shallow copy
 * of an extant #ObitImageMosaicWB.  
 *
 * A copy of a pointer to an ObitImageMosaicWB should always be made using the
 * ObitImageMosaicWBRef function which updates the reference count in the object.
 * Then whenever freeing an ObitImageMosaicWB or changing a pointer, the function
 * ObitImageMosaicWBUnref will decrement the reference count and destroy the 
 * object when the reference count hits 0.
 *
 */

/*--------------Class definitions-------------------------------------*/
/** ObitImageMosaicWB Class structure. */
typedef struct {
#include "ObitImageMosaicWBDef.h"   /* this class definition */
} ObitImageMosaicWB;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitImageMosaicWB
 * returns a ObitImageMosaicWB*.
 * in = object to unreference
 */
#define ObitImageMosaicWBUnref(in) ObitUnref ( in)

/** 
 * Macro to reference (update reference count) an ObitImageMosaicWB.
 * returns a ObitImageMosaicWB*.
 * in = object to reference
 */
#define ObitImageMosaicWBRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitImageMosaicWBIsA(in) ObitIsA (in, ObitImageMosaicWBGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitImageMosaicWBClassInit (void);

/** Public: Constructor. */
ObitImageMosaicWB* newObitImageMosaicWB (gchar* name, olong number);

/** Public: Create ImageMosaicWB object from description in an ObitInfoList */
ObitImageMosaicWB* ObitImageMosaicWBFromInfo (gchar *prefix, ObitInfoList *inList, 
					      ObitErr *err);
typedef ObitImageMosaicWB* (*ObitImageMosaicWBFromInfoFP) (gchar *prefix, 
							   ObitInfoList *inList, 
							   ObitErr *err);
/** Public: ClassInfo pointer */
gconstpointer ObitImageMosaicWBGetClass (void);

/** Public: Copy (shallow) constructor. */
ObitImageMosaic* 
ObitImageMosaicWBCopy  (ObitImageMosaic *in, ObitImageMosaic *out, 
			ObitErr *err);

/** Public: Zap specified image. */
void 
ObitImageMosaicWBZapImage  (ObitImageMosaic *in, olong number,
			    ObitErr *err);

/** Public: Set underlying files */
void ObitImageMosaicWBSetFiles  (ObitImageMosaic *in, gboolean doBeam, ObitErr *err);

/** Public: Create Mosaic from uv data */
ObitImageMosaicWB *ObitImageMosaicWBCreate (gchar *name, olong order, ObitUV *uvData, 
					    ObitErr *err);

/** Public: Define parameters of images */
void ObitImageMosaicWBDefine (ObitImageMosaic *in, ObitUV *uvData, gboolean doBeam,
			      ObitErr *err);

/** Public: Flatten tiles onto full field image */
void ObitImageMosaicWBFlatten (ObitImageMosaic *in, ObitErr *err);

/** Public: Reimaging needed to center strong source on pixel? */
gboolean ObitImageMosaicWBReimage (ObitImageMosaicWB *mosaic, ObitErr* err);

/** Public: Add field to mosaic */
void ObitImageMosaicWBAddField (ObitImageMosaic *in, ObitUV *uvData, 
				olong nx, olong ny, olong nplane, 
				ofloat RAShift, ofloat DecShift, 
				gboolean isAuto, ObitErr *err);

/** Public:  Generate a mosaic for peeling */
ObitImageMosaicWB* ObitImageMosaicWBMaxField (ObitImageMosaic *mosaic, 
					      ofloat MinFlux, olong *ignore, olong *field,
					      ObitErr* err); 

/** Public: Extract information about underlying structures to ObitInfoList */
void ObitImageMosaicWBGetInfo (ObitImageMosaic *in, gchar *prefix, ObitInfoList *outList, 
			       ObitErr *err);

/** Public: Concatenate Image CC tables onto the FullField Image */
void ObitImageMosaicWBCopyCC (ObitImageMosaic *in, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitImageMosaicWBClassDef.h"
} ObitImageMosaicWBClassInfo; 

#endif /* OBITIMAGEMOSAICWB_H */ 
