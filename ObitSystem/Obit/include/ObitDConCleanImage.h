/* $Id: ObitDConCleanImage.h,v 1.5 2007/08/31 17:24:48 bcotton Exp $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                          */
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
#ifndef OBITDCONCLEANIMAGE_H 
#define OBITDCONCLEANIMAGE_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitDConClean.h"
#include "ObitImageMosaic.h"
#include "ObitDConCleanWindow.h"
#include "ObitDConCleanBmHist.h"
#include "ObitDConCleanPxHist.h"
#include "ObitDConCleanPxList.h"
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanImage.h
 * ObitDConCleanImage Image based (BGC-like) CLEAN class.
 *
 * This class is derived from the #ObitDConClean class.
 * 
 * Only a Single field can be CLEANed and only a quarter the area of it.
 *
 * \section ObitDConCleanImageaccess Creators and Destructors
 * An ObitDConCleanImage will usually be created using ObitDConCleanImageCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitDConCleanImage should always be made using the
 * #ObitDConCleanImageRef function which updates the reference count in the object.
 * Then whenever freeing an ObitDConCleanImage or changing a pointer, the function
 * #ObitDConCleanImageUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 * \section ObitDConCleanImagecontrol CLEAN control information
 * The control parameters for the CLEAN are read from the ObitInfoList member
 * when the Deconvolve function is called:
 * The file etc. info should have been stored in the ObitInfoList:
 * \li "Niter"   OBIT_long scalar   = Maximum number of CLEAN iterations
 * \li "maxPixel" OBIT_long scalar  = Maximum number of residuals [def 20000]
 * \li "minPatch" OBIT_long scalar  = Minimum beam patch in pixels [def 100]
 * \li "BMAJ"    OBIT_float scalar = Restoring beam major axis (deg)
 * \li "BMIN"    OBIT_float scalar = Restoring beam minor axis (deg)
 * \li "BPA"     OBIT_float scalar = Restoring beam position angle (deg)
 * \li "Beam"    OBIT_float array[3]= (BMAJ, BMIN, BPA) alternate form
 * \li "Gain"    OBIT_float scalar = CLEAN loop gain
 * \li "minFlux" OBIT_float scalar = Minimun flux density (Jy)
 * \li "Factor"  OBIT_float array  = CLEAN depth factor
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 *                                   [def 1,1,1,1,1]
 * \li "CLEANBox" OBIT_long [4,?]   = Array of Clean boxes for field 1 
 *                                   Any entries with first element=0 are ignored.
 * \li "prtLv"   OBIT_int            message level  [def 2]
 *                       0=none, 1=summary, 2=normal, higher numbers for diagnostics
 * \li "CCVer"   OBIT_long          = CC table version number
 */

/*--------------Class definitions-------------------------------------*/
/** ObitDConCleanImage Class structure. */
typedef struct {
#include "ObitDConCleanImageDef.h"   /* this class definition */
} ObitDConCleanImage;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitDConCleanImage
 * returns a ObitDConCleanImage*.
 * in = object to unreference
 */
#define ObitDConCleanImageUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitDConCleanImage.
 * returns a ObitDConCleanImage*.
 * in = object to reference
 */
#define ObitDConCleanImageRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitDConCleanImageIsA(in) ObitIsA (in, ObitDConCleanImageGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitDConCleanImageClassInit (void);

/** Public: Default Constructor. */
ObitDConCleanImage* newObitDConCleanImage (gchar* name);

/** Public: Create/initialize ObitDConCleanImage structures */
ObitDConCleanImage* ObitDConCleanImageCreate (gchar* name, ObitImageMosaic *mosaic, 
					      ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitDConCleanImageGetClass (void);

/** Public: Copy (deep) constructor. */
ObitDConCleanImage* ObitDConCleanImageCopy  (ObitDConCleanImage *in, 
					     ObitDConCleanImage *out, ObitErr *err);

/** Public: Copy structure. */
void ObitDConCleanImageClone (ObitDConCleanImage *in, 
			      ObitDConCleanImage *out, ObitErr *err);

/** Public: Do deconvolution. */
void ObitDConCleanImageDeconvolve (ObitDCon *in, ObitErr *err);

/** Public: Get parameters. */
void  ObitDConCleanImageGetParms (ObitDCon *in, ObitErr *err);

/** Public: Subtract components and generate new residual image(s). */
void ObitDConCleanImageSub(ObitDConCleanImage *in, ObitErr *err);


/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitDConCleanImageClassDef.h"
} ObitDConCleanImageClassInfo; 

#endif /* OBITDCONCLEANIMAGE_H */ 
