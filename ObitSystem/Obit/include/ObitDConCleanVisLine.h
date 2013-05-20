/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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
#ifndef OBITDCONCLEANVISLINE_H 
#define OBITDCONCLEANVISLINE_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitDConCleanPxList.h"
#include "ObitDConClean.h"
#include "ObitImageMosaic.h"
#include "ObitDConCleanWindow.h"
#include "ObitDConCleanBmHist.h"
#include "ObitDConCleanPxHist.h"
#include "ObitUVImager.h"
#include "ObitSkyModel.h"
#include "ObitDisplay.h"
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanVisLine.h
 * ObitDConCleanVisLine Visibility-based (Cotton-Schwab) 
 * Spectral line CLEAN class.
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
 * Parallel processing of spectral line datasets are supported.
 *
 * \section ObitDConCleanVisLineaccess Creators and Destructors
 * An ObitDConCleanVisLine will usually be created using ObitDConCleanVisLineCreate 
 * or ObitDConCleanVisLineCreate2 which allow
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitDConCleanVisLine should always be made using the
 * #ObitDConCleanVisLineRef function which updates the reference count in the object.
 * Then whenever freeing an ObitDConCleanVisLine or changing a pointer, the function
 * #ObitDConCleanVisLineUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 * \section ObitDConCleanVisLinecontrol CLEAN control information
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
/** ObitDConCleanVisLine Class structure. */
typedef struct {
#include "ObitDConCleanVisLineDef.h"   /* this class definition */
} ObitDConCleanVisLine;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitDConCleanVisLine
 * returns a ObitDConCleanVisLine*.
 * in = object to unreference
 */
#define ObitDConCleanVisLineUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitDConCleanVisLine.
 * returns a ObitDConCleanVisLine*.
 * in = object to reference
 */
#define ObitDConCleanVisLineRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitDConCleanVisLineIsA(in) ObitIsA (in, ObitDConCleanVisLineGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitDConCleanVisLineClassInit (void);

/** Public: Default Constructor. */
ObitDConCleanVisLine* newObitDConCleanVisLine (gchar* name);

/** Public: Create/initialize ObitDConCleanVisLine structures */
ObitDConCleanVisLine* 
ObitDConCleanVisLineCreate (gchar* name, olong nPar, olong nAvg, ObitUV *uvdata, 
			    ObitErr *err);

/** Public: Create/initialize ObitDConCleanVisLine structures from
    optional components */
ObitDConCleanVisLine* 
ObitDConCleanVisLineCreate2 (gchar* name, ObitUV *uvdata, 
			     ObitUVImager *imager, ObitSkyModel *skyModel, 
			     ObitErr *err);

/** Public: Update channel selection */
void ObitDConCleanVisLineUpdate (ObitDConCleanVisLine *in, olong nPar, 
				 olong BIF, olong BChan);

/** Public: ClassInfo pointer */
gconstpointer ObitDConCleanVisLineGetClass (void);

/** Public: Copy (deep) constructor. */
ObitDConCleanVisLine* ObitDConCleanVisLineCopy (ObitDConCleanVisLine *in, 
						ObitDConCleanVisLine *out, ObitErr *err);

/** Public: Copy structure. */
void ObitDConCleanVisLineClone (ObitDConCleanVisLine *in, 
			    ObitDConCleanVisLine *out, ObitErr *err);

/** Public: Do deconvolution. */
void ObitDConCleanVisLineDeconvolve (ObitDCon *in, ObitErr *err);

/** Public: Get parameters. */
void  ObitDConCleanVisLineGetParms (ObitDCon *in, ObitErr *err);

/** Public: Subtract components from UV data */
void ObitDConCleanVisLineSub(ObitDConClean *in, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitDConCleanVisLineClassDef.h"
} ObitDConCleanVisLineClassInfo; 

#endif /* OBITDCONCLEANVISLINE_H */ 
