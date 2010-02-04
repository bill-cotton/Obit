/* $Id: ObitDConCleanVisWB.h 128 2009-09-23 14:48:29Z bill.cotton $   */
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
#ifndef OBITDCONCLEANVISWB_H 
#define OBITDCONCLEANVISWB_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitDConCleanPxList.h"
#include "ObitDConClean.h"
#include "ObitDConCleanVis.h"
#include "ObitImageMosaic.h"
#include "ObitImageMosaicWB.h"
#include "ObitDConCleanWindow.h"
#include "ObitDConCleanBmHist.h"
#include "ObitDConCleanPxHist.h"
#include "ObitUVImager.h"
#include "ObitUVImagerWB.h"
#include "ObitSkyModel.h"
#include "ObitDisplay.h"
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanVisWB.h
 * ObitDConCleanVisWB Visibility-based (Cotton-Schwab) CLEAN class 
 * for wideband imaging in the style of Sault & Wieringa 1994 A&AS 108, 585
 * Fit S = S_0 exp(alpha log(nu/nu_0) + beta log(nu/nu_0)*2)
 * using approximation:
 * S = S0(1+ alpha*log(nu/nu_0) + beta * log(nu/nu_0)^2)
 *
 * This class is derived from the #ObitDConCleanVis class.
 *
 * autoWindow feature will automatically set CLEAN windows inside 
 * a predefined outer window.  Each cycle the residuals inside the outer 
 * window are searched to the maximum value; if the peak is outside the 
 * inner window and > 3 sigma, a new round box of radius 3 pixels is added 
 * to the window.  Cleaning in each cycle will stop when the peak residual 
 * drops to the level of the highest value outside the CLEAN window.
 * In autoWindow mode, the field to be processed next is determined using 
 * the statistics from the outer window and for normal mode, the inner window.
 * This should result in the brightest emission being cleaned next 
 * and a box added on it if necessary.
 *
 * \section ObitDConCleanVisWBaccess Creators and Destructors
 * An ObitDConCleanVisWB will usually be created using ObitDConCleanVisWBCreate 
 * or ObitDConCleanVisWBCreate2 which allow
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitDConCleanVisWB should always be made using the
 * #ObitDConCleanVisWBRef function which updates the reference count in the object.
 * Then whenever freeing an ObitDConCleanVisWB or changing a pointer, the function
 * #ObitDConCleanVisWBUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 * \section ObitDConCleanVisWBcontrol CLEAN control information
 * The control parameters for the CLEAN are read from the ObitInfoList member
 * when the Deconvolve function is called:
 * The file etc. info should have been stored in the ObitInfoList:
 * \li "Niter"   OBIT_long scalar   = Maximum number of CLEAN iterations
 * \li "maxPixel" OBIT_long scalar  = Maximum number of residuals [def 20000]
 * \li "minPatch" OBIT_long scalar  = Minimum beam patch in pixels [def 100]
 * \li "BMAJ"    OBIT_float scalar = Restoring beam major axis (deg)
 * \li "BMIN"    OBIT_float scalar = Restoring beam minor axis (deg)
 * \li "BPA"     OBIT_float scalar = Restoring beam position angle (deg) (",", deg)
 * \li "Beam"    OBIT_float array[3]= (BMAJ, BMIN, BPA) alternate form
 * \li "Gain"    OBIT_float scalar = CLEAN loop gain per field
 *                                   If one given, used for all
 * \li "minFlux" OBIT_float scalar = Minimum flux density (Jy) per field
 *                                   If one given, used for all
 * \li "Factor"  OBIT_float array  = CLEAN depth factor per field
 *                                   If one given, used for all
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 * \li "CCVer"   OBIT_long          = CC table version number
 * \li "CLEANBox"OBIT_long [4,?]    = Array of Clean boxes for field 1
 *                                   Any entries with first element=0 are ignored.
 * \li "prtLv"    OBIT_int           message level  [def 2]
 *                       0=none, 1=summary, 2=normal, higher numbers for diagnostics
 * \li "autoWindow" OBIT_boolean scalar = True if autoWindow feature wanted.
 * \li "Mode"      OBIT_long scalar = Model mode (ObitSkyModelMode) [def OBIT_SkyModel_Fastest]
 * \li "doRestore" OBIT_bool       = Restore image when done? [def TRUE]
 * \li "doXRestore" OBIT_bool      = Cross restore images when done? [def doRestore]
 * \li "doFlatten" OBIT_bool       = Flatten image when done? [def TRUE]
 * \li "doWeight"  OBIT_bool       = Weight UV data before imaging? [def TRUE]
 * \li "doBeam"    OBIT_bool       = Need to (re)make beam? [def TRUE]
 * \li "doRecenter" OBIT_bool      = Allow recentering autoCenter fields? [def TRUE]
 * \li "reuseFlux" OBIT_float      = Level of Components in initial 
 *                                   SkyModel to reuse, <0 -> none [def none]
 * \li "autoCen" OBIT_float scalar = Threshold leven for autocenter
 *               If peak exceeds this value the minFlux reset to 0.1*autoCen
 * \li "dispURL"   OBIT_string     = URL of display server
 */

/*--------------Class definitions-------------------------------------*/
/** ObitDConCleanVisWB Class structure. */
typedef struct {
#include "ObitDConCleanVisWBDef.h"   /* this class definition */
} ObitDConCleanVisWB;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitDConCleanVisWB
 * returns a ObitDConCleanVisWB*.
 * in = object to unreference
 */
#define ObitDConCleanVisWBUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitDConCleanVisWB.
 * returns a ObitDConCleanVisWB*.
 * in = object to reference
 */
#define ObitDConCleanVisWBRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitDConCleanVisWBIsA(in) ObitIsA (in, ObitDConCleanVisWBGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitDConCleanVisWBClassInit (void);

/** Public: Default Constructor. */
ObitDConCleanVisWB* newObitDConCleanVisWB (gchar* name);

/** Public: Create/initialize ObitDConCleanVisWB structures */
ObitDConCleanVisWB* ObitDConCleanVisWBCreate (gchar* name, ObitUV *uvdata, 
					      olong order, ObitErr *err);

/** Public: Create/initialize ObitDConCleanVisWB structures from
    optional components */
ObitDConCleanVis* 
ObitDConCleanVisWBCreate2 (gchar* name, ObitUV *uvdata, 
			   ObitUVImager *imager, ObitSkyModel *skyModel, 
			   olong order, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitDConCleanVisWBGetClass (void);

/** Public: Copy (deep) constructor. */
ObitDConCleanVis* ObitDConCleanVisWBCopy  (ObitDConCleanVis *in, 
					   ObitDConCleanVis *out, ObitErr *err);

/** Public: Copy structure. */
void ObitDConCleanVisWBClone (ObitDConCleanVis *in, 
			      ObitDConCleanVis *out, ObitErr *err);

/** Public: Do deconvolution. */
void ObitDConCleanVisWBDeconvolve (ObitDCon *in, ObitErr *err);

/** Public: Get parameters. */
void  ObitDConCleanVisWBGetParms (ObitDCon *in, ObitErr *err);

/** Public: Set Default CLEAN windows */
void  ObitDConCleanVisWBDefWindow (ObitDConClean *in, ObitErr *err);

/** Public: Subtract components from UV data */
void ObitDConCleanVisWBSub(ObitDConClean *in, ObitErr *err);

/** Public: Determine quality measure for field */
ofloat ObitDConCleanVisWBQuality(ObitDConCleanVis *in, olong field, 
				 ObitErr *err);

/** Public: Reimaging needed to center strong source on pixel? */
gboolean ObitDConCleanVisWBReimage (ObitDConCleanVisWB *in, ObitUV* uvdata, 
				  ObitErr* err);

/** Public: Resize to add a new field */
void ObitDConCleanVisWBAddField (ObitDConCleanVis *in, ObitUV* uvdata, 
				 ObitErr* err);

/** Public: Recenter autoCenter images */
gboolean ObitDConCleanVisWBRecenter (ObitDConCleanVis *in, ObitUV* uvdata, 
				     ObitErr* err);

/** Public: Filter weak, isolated components */
gboolean ObitDConCleanVisWBFilter (ObitDConCleanVis *in, ofloat filter[2], 
				   ObitErr* err);

/** Public: Automatically add window. */
gboolean ObitDConCleanVisWBAutoWindow(ObitDConClean *in, olong *fields, 
				      ObitFArray **pixarray, ObitErr *err);

/** Public:  Get pixel statistics. */
gboolean ObitDConCleanVisWBPixelStats(ObitDConClean *in, ObitFArray **pixarray,
				      ObitErr *err);

/** Public: Select CLEAN components? */
gboolean ObitDConCleanVisWBSelect(ObitDConClean *inn, ObitFArray **pixarray, 
				  ObitErr *err);

/** Public:  Restore subtracted components. */
void ObitDConCleanVisWBRestore(ObitDConClean *in, ObitErr *err);

/** Public:  Restore subtracted components from other fields. */
void ObitDConCleanVisWBXRestore(ObitDConClean *in, ObitErr *err);

/** Public:  Normalize spectral images. */
void ObitDConCleanVisWBSpecNorm(ObitDConCleanVisWB *in, ObitErr *err);

/** Public: Flatten multiple facets to one. */
void ObitDConCleanVisWBFlatten(ObitDConClean *in, ObitErr *err);

/** Public:  Decompose images into spectral components. */
void ObitDConCleanVisWBDecomp(ObitDConCleanVisWB *in, ObitErr *err);
typedef void (*ObitDConCleanVisWBDecompFP) (ObitDConCleanVisWB *in, 
					    ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitDConCleanVisWBClassDef.h"
} ObitDConCleanVisWBClassInfo; 

#endif /* OBITDCONCLEANVISWB_H */ 
