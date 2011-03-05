/* $Id$       */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2011                                          */
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
#ifndef OBITDCONCLEAN_H 
#define OBITDCONCLEAN_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitDCon.h"
#include "ObitImageMosaic.h"
#include "ObitDConCleanWindow.h"
#include "ObitDConCleanBmHist.h"
#include "ObitDConCleanPxHist.h"
#include "ObitDConCleanPxList.h"
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConClean.h
 * ObitDConClean virtual CLEAN base class.
 *
 * This class is derived from the #ObitDCon class.
 *
 * Actual deconvolution classes are derived from this class
 * autoWindow feature will automatically set CLEAN windows inside 
 * a predefined outer window.  Each cycle the residuals inside the outer 
 * window are searched to the maximum value; if the peak is outside the 
 * inner window and > 3 sigma, a new round box of radius 3 pixels is added 
 * to the window.  Cleaning in each cycle will stop when the peak residual 
 * drops to the level of the highest value outside the CLEAN window.
 * 
 * \section ObitDConCleanaccess Creators and Destructors
 * An ObitDConClean will usually be created using ObitDConCleanCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitDConClean should always be made using the
 * #ObitDConCleanRef function which updates the reference count in the object.
 * Then whenever freeing an ObitDConClean or changing a pointer, the function
 * #ObitDConCleanUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 * \section ObitDConCleancontrol CLEAN control information
 * The control parameters for the CLEAN are read from the ObitInfoList member
 * when the Deconvolve function is called:
 * The file etc. info should have been stored in the ObitInfoList:
 * \li "Niter"   OBIT_long scalar   = Maximum number of CLEAN iterations
 * \li "maxPixel" OBIT_long scalar  = Maximum number of residuals [def 20000]
 * \li "minPatch" OBIT_long scalar  = Minimum beam patch in pixels [def 100]
 * \li "BMAJ"    OBIT_float scalar = Restoring beam major axis (deg)
 * \li "BMIN"    OBIT_float scalar = Restoring beam minor axis (deg)
 * \li "BPA"     OBIT_float scalar = Restoring beam position angle (deg)
 * \li "Beam"    OBIT_float array[3]= (BMAJ, BMIN, BPA) alternate form (",", deg)
 * \li "Gain"    OBIT_float array  = CLEAN loop gain per field
 *                                   If only one given it is used for all.
 * \li "minFlux" OBIT_float array  = Minimum flux density (Jy)  per field
 *                                   If only one given it is used for all.
 * \li "Factor"  OBIT_float array  = CLEAN depth factor per field
 *                                   If only one given it is used for all.
 * \li "CCVer"   OBIT_long          = CC table version number
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 *                                   def (1,1,1,1,1)
 * \li "CLEANBox" OBIT_long [4,?]   = Array of Clean boxes for field 1
 *                                   Any entries with first element=0 are ignored.
 * \li "prtLv"    OBIT_int           message level  [def 2]
 *                       0=none, 1=summary, 2=normal, higher numbers for diagnostics
 * \li "autoWindow" OBIT_boolean scalar = True if autoWindow feature wanted.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitDConClean Class structure. */
typedef struct {
#include "ObitDConCleanDef.h"   /* this class definition */
} ObitDConClean;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitDConClean
 * returns a ObitDConClean*.
 * in = object to unreference
 */
#define ObitDConCleanUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitDConClean.
 * returns a ObitDConClean*.
 * in = object to reference
 */
#define ObitDConCleanRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitDConCleanIsA(in) ObitIsA (in, ObitDConCleanGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitDConCleanClassInit (void);

/** Public: Default Constructor. */
ObitDConClean* newObitDConClean (gchar* name);

/** Public: Create/initialize ObitDConClean structures */
ObitDConClean* ObitDConCleanCreate (gchar* name, ObitImageMosaic *mosaic, 
			  ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitDConCleanGetClass (void);

/** Public: Copy (deep) constructor. */
ObitDConClean* ObitDConCleanCopy  (ObitDConClean *in, ObitDConClean *out, ObitErr *err);

/** Public: Copy structure. */
void ObitDConCleanClone (ObitDConClean *in, ObitDConClean *out, ObitErr *err);

/** Public: Do deconvolution. */
void ObitDConCleanDeconvolve (ObitDCon *in, ObitErr *err);

/** Public: Get parameters. */
void  ObitDConCleanGetParms (ObitDCon *in, ObitErr *err);

/** Public: Set Default CLEAN windows. */
void ObitDConCleanDefWindow(ObitDConClean *in, ObitErr *err);
typedef void (*ObitDConCleanDefWindowFP) (ObitDConClean *in, ObitErr *err);

/** Public:  Prepare for minor cycle. */
gboolean ObitDConCleanPixelStats(ObitDConClean *in, ObitFArray **pixarray,
				 ObitErr *err);
typedef gboolean (*ObitDConCleanPixelStatsFP) (ObitDConClean *in, ObitFArray **pixarray,
					       ObitErr *err);

/** Public:  Determine image statistics. */
void ObitDConCleanImageStats(ObitDConClean *in, olong field, gboolean doBeam, 
			     ObitErr *err);
typedef void (*ObitDConCleanImageStatsFP) (ObitDConClean *in, olong field, 
					   gboolean doBeam, ObitErr *err);
 
/** Public:Select components to be subtracted . */
gboolean ObitDConCleanSelect(ObitDConClean *in, ObitFArray **pixarray, ObitErr *err);
typedef gboolean (*ObitDConCleanSelectFP) (ObitDConClean *in, ObitFArray **pixarray,
					   ObitErr *err);

/** Public: Subtract components and generate new residual image(s). */
void ObitDConCleanSub(ObitDConClean *in, ObitErr *err);
typedef void (*ObitDConCleanSubFP) (ObitDConClean *in, ObitErr *err);

/** Public:  Restore subtracted components. */
void ObitDConCleanRestore(ObitDConClean *in, ObitErr *err);
typedef void (*ObitDConCleanRestoreFP) (ObitDConClean *in, ObitErr *err);

/** Public:  Restore subtracted components from other fields. */
void ObitDConCleanXRestore(ObitDConClean *in, ObitErr *err);
typedef void (*ObitDConCleanXRestoreFP) (ObitDConClean *in, ObitErr *err);

/** Public: Flatten multiple facets to one. */
void ObitDConCleanFlatten(ObitDConClean *in, ObitErr *err);
typedef void (*ObitDConCleanFlattenFP) (ObitDConClean *in, ObitErr *err);

/** Public: Automatically add window. */
gboolean ObitDConCleanAutoWindow(ObitDConClean *in, olong *fields, 
				 ObitFArray **pixarray, ObitErr *err);
typedef gboolean (*ObitDConCleanAutoWindowFP) (ObitDConClean *in, olong *fields, 
					       ObitFArray **pixarray, ObitErr *err);

/** Public: Get cross restoring beam parameters */
ofloat ObitDConCleanGetXRestoreBeam(ObitImageDesc *imDesc1, ObitImageDesc *imDesc2, 
				    ofloat *gparm, ofloat *bmaj, ofloat *bmin, ofloat *bpa);

/* Private routines */
/** Private: Read Beam patches. */
void ReadBP (ObitDConClean* in, ObitErr *err);
typedef void (*ReadBPFP) (ObitDConClean* in, ObitErr *err);


/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitDConCleanClassDef.h"
} ObitDConCleanClassInfo; 

#endif /* OBITDCONCLEAN_H */ 
