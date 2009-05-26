/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005,2009                                          */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITDCONCLEANOTF_H 
#define OBITDCONCLEANOTF_H 

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
 * \file ObitDConCleanOTF.h
 * ObitDConCleanOTF Image based (BGC-like) CLEAN class for OTF (single dish data).
 * This CLEAN is appropriate for Single dish deconvolution where the 
 * instrumental PSF is of very limited support.
 * The dirty image is used with the dirty beam to decompose selected pixels into a
 * set of delta functions stores in an AIPS CC table.
 * Windowing is supported and the dirty beam should have the same pixel spacing 
 * as the dirty beam but need not be the same size.
 *
 * \section ObitDConCleanOTFaccess Creators and Destructors
 * An ObitDConCleanOTF will usually be created using ObitDConCleanOTFCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitDConCleanOTF should always be made using the
 * #ObitDConCleanOTFRef function which updates the reference count in the object.
 * Then whenever freeing an ObitDConCleanOTF or changing a pointer, the function
 * #ObitDConCleanOTFUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 * \section ObitDConCleanOTFcontrol CLEAN control information
 * The control parameters for the CLEAN are read from the ObitInfoList member
 * when the Deconvolve function is called:
 * The file etc. info should have been stored in the ObitInfoList:
 * \li "Niter"   OBIT_long scalar   = Maximum number of CLEAN iterations
 * \li "Patch"   OBIT_long scalar   = Size of beam patch in pixels [def all]
 * \li "BeamSize"OBIT_float scalar = Restoring beam FWHM (deg)
 * \li "Gain"    OBIT_float scalar = CLEAN loop gain
 * \li "minFlux" OBIT_float scalar = Minimun flux density (Jy)
 * \li "Factor"  OBIT_float array  = CLEAN depth factor
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 * \li "CCVer"   OBIT_long          = CC table version number
 * \li "noResid" OBIT_bool scalar  = If defined and TRUE, set residuals to Zero in restore
 * \li doScale   OBIT_bool scalar If defined and TRUE, Scale residuals by 
 *               ratio of beam areas, def TRUE
 * \li "doScaleCC" OBIT_bool scalar = Scale CCs by ratio of dirty to CLEAN beam areas?
 */

/*--------------Class definitions-------------------------------------*/
/** ObitDConCleanOTF Class structure. */
typedef struct {
#include "ObitDConCleanOTFDef.h"   /* this class definition */
} ObitDConCleanOTF;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitDConCleanOTF
 * returns a ObitDConCleanOTF*.
 * in = object to unreference
 */
#define ObitDConCleanOTFUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitDConCleanOTF.
 * returns a ObitDConCleanOTF*.
 * in = object to reference
 */
#define ObitDConCleanOTFRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitDConCleanOTFIsA(in) ObitIsA (in, ObitDConCleanOTFGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitDConCleanOTFClassInit (void);

/** Public: Default Constructor. */
ObitDConCleanOTF* newObitDConCleanOTF (gchar* name);

/** Public: Create/initialize ObitDConCleanOTF structures */
ObitDConCleanOTF* 
ObitDConCleanOTFCreate (gchar* name, ObitImage *dirty, ObitImage *beam,  
			ObitImage *clean, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitDConCleanOTFGetClass (void);

/** Public: Copy (deep) constructor. */
ObitDConCleanOTF* ObitDConCleanOTFCopy  (ObitDConCleanOTF *in, 
					 ObitDConCleanOTF *out, ObitErr *err);

/** Public: Copy structure. */
void ObitDConCleanOTFClone (ObitDConCleanOTF *in, 
			    ObitDConCleanOTF *out, ObitErr *err);

/** Public: Do deconvolution. */
void ObitDConCleanOTFDeconvolve (ObitDCon *in, ObitErr *err);

/** Public: Get parameters. */
void  ObitDConCleanOTFGetParms (ObitDCon *in, ObitErr *err);

/** Public: Generate residual image */
void ObitDConCleanOTFSub(ObitDConClean *in, ObitErr *err);

/** Public:Select components  . */
gboolean ObitDConCleanOTFSelect(ObitDConClean *in, ObitFArray *pixarray,
				ObitErr *err);

/** Public: Subtract components and generate new residual image(s). */
void ObitDConCleanOTFSub(ObitDConClean *in, ObitErr *err);

/** Public:  Restore subtracted components. */
void ObitDConCleanOTFRestore(ObitDConClean *in, ObitErr *err);

/** Public:  Scale CCs by ratio of dirty to CLEAN beam areas. */
void ObitDConCleanOTFScaleCC(ObitDConCleanOTF *in, ObitErr *err);
typedef void (*ObitDConCleanOTFScaleCCFP) (ObitDConCleanOTF *in, ObitErr *err);

/** Public:  Prepare for minor cycle. */
void ObitDConCleanOTFPixelStats(ObitDConClean *in, ObitFArray *pixarray, 
				ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitDConCleanOTFClassDef.h"
} ObitDConCleanOTFClassInfo; 

#endif /* OBITDCONCLEANOTF_H */ 
