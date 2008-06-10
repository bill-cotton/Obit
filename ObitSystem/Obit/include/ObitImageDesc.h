/* $Id: ObitImageDesc.h,v 1.7 2007/08/31 17:24:48 bcotton Exp $   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITIMAGEDESC_H 
#define OBITIMAGEDESC_H 
#include "Obit.h"
#include "ObitErr.h"
#include "ObitInfoList.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageDesc.h
 * ObitImageDesc Obit Image descriptor class definition.
 *
 * This class is derived from the #Obit class.
 *
 * This contains information about the observations and the size and 
 * coordinates in the image.
 * Also included are the current location of the image in an #ObitImage
 * image buffer and the specified subimaging parameters.
 *
 * \section ObitImageDescUsage Usage
 * Instances can be obtained using the #newObitImageDesc constructor
 * the #ObitImageDescCopy copy constructor or a pointer duplicated using 
 * the #ObitImageDescRef function.
 * When an instance is no longer needed, use the #ObitImageDescUnref macro
 * to release it.
 */

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitImageDesc
 * returns a ObitImageDesc* (NULL).
 * \li in = object to unreference.
 */
#define ObitImageDescUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitImageDesc.
 * returns a ObitImageDesc*.
 * in = object to reference
 */
#define ObitImageDescRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitImageDescIsA(in) ObitIsA (in, ObitImageDescGetClass())

/** Maximum number of dimensions */
#define IM_MAXDIM 7       /* maximum array dimension */
/** Maximum length of descriptor string value */
#define IMLEN_VALUE 41
/** Maximum length of descriptor keyword  */
#define IMLEN_KEYWORD 21

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitCoordType
 * enum for coordinate system types (Equatorial, Galactic, Ecliptic)
 */
enum obitCoordType {
  /** Equatorial (RA, Dec) */
  OBIT_Equatorial, 
  /** Galactic (GLONG, GLAT) */
  OBIT_Galactic,  
  /** Ecliptic (ELONG, ELAT) */
  OBIT_Ecliptic
}; /* end enum obitCoordType */
/** typedef for enum for coordinate system types. */
typedef enum obitCoordType ObitCoordType;
/*--------------Class definitions-------------------------------------*/
/**
 * ObitImageDesc Class structure.
 *
 * This class contains descriptions of astronomical images.
 */  
typedef struct {
#include "ObitImageDescDef.h"  /* actual definitions */
} ObitImageDesc;

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitImageDescClassInit (void);

/** Public: Constructor. */
ObitImageDesc* newObitImageDesc (gchar *name);

/** Public: Return class pointer. */
gconstpointer ObitImageDescGetClass (void);

/** Public: Copy ImageDesc */
ObitImageDesc* 
ObitImageDescCopy (ObitImageDesc* in, 
		   ObitImageDesc* out, ObitErr *err);

/** Public: Copy descriptive (nonstructural) information. */
void ObitImageDescCopyDesc (ObitImageDesc* in, ObitImageDesc* out,
			 ObitErr *err);

/** Public: Default ImageDesc */
ObitImageDesc* ObitImageDescDefault (gchar *name);

/** Public: Transform pixel in one image to another */
void ObitImageDescCo (ObitImageDesc* in, ObitImageDesc* out,
			 ObitErr *err);

/** Public: Index for easier access */
void ObitImageDescIndex (ObitImageDesc* in);

/** Public: Convert pixels on one image to those in another */
gboolean ObitImageDescCvtPixel(ObitImageDesc* in, ObitImageDesc* out, 
			       ofloat *inPixel, ofloat *outPixel, ObitErr *err);

/** Public: Convert pixels on one image to those in another with 
    Zernike corrections*/
gboolean 
ObitImageDescCvtZern(ObitImageDesc* in, ObitImageDesc* out, 
		     olong nZern, ofloat *ZCoef,
		     ofloat *inPixel, ofloat *outPixel, ObitErr *err);

/** Public: Determine the location of a pixel in an image */
void ObitImageDescGetPos(ObitImageDesc* in, ofloat *inPixel, 
			 odouble *pos, ObitErr *err);

/** Public: Determine the pixel of a location in an image */
void ObitImageDescGetPixel(ObitImageDesc* in, odouble *pos, 
			   ofloat *outPixel, ObitErr *err);

/**  Public: Is there overlap in two images */
gboolean ObitImageDescOverlap(ObitImageDesc *in1, ObitImageDesc *in2, ObitErr *err);

/**  Public: Tell rotation angle of image */
ofloat ObitImageDescRotate(ObitImageDesc *in);

/**  Public: Tell Pointing position */
void ObitImageDescGetPoint(ObitImageDesc *in, odouble *RAPnt, odouble *DecPnt);

/** Public:  Determine the angular distance to the antenna pointing position */
ofloat ObitImageDescAngle (ObitImageDesc *in, ofloat y, ofloat x);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitImageDescClassDef.h" /* Actual definition */
} ObitImageDescClassInfo; 

#endif /* OBITIMAGEDESC_H */ 

