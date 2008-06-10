/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2008                                          */
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
#ifndef OBITDCONCLEANOTFREC_H 
#define OBITDCONCLEANOTFREC_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitOTF.h"
#include "ObitDConClean.h"
#include "ObitDisplay.h"
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanOTFRec.h
 *
 * This class derived from the #ObitDConClean class.
 *
 * ObitDConCleanOTF Record based CLEAN class for ObitOTF (single dish data).
 * The technique is to occasionally subtract the CLEAN components from
 * the OTF data and reimage.  This is most useful for beamswitched data.
 * The dirty image is used with the dirty beam to decompose selected pixels into a
 * set of delta functions stored in an AIPS CC table.
 * Windowing is supported and the dirty beam should have the same pixel spacing 
 * as the dirty beam but need not be the same size.
 *
 * \section ObitDConCleanOTFaccess Creators and Destructors
 * An ObitDConCleanOTF will usually be created using ObitDConCleanOTFCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitDConCleanOTFRec should always be made using the
 * #ObitDConCleanOTFRecRef function which updates the reference count in the object.
 * Then whenever freeing an ObitDConCleanOTFRec or changing a pointer, the function
 * #ObitDConCleanOTFRecUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 * \section ObitDConCleanOTFcontrol CLEAN control information
 * The control parameters for the CLEAN are read from the ObitInfoList member
 * when the Deconvolve function is called:
 * The file etc. info should have been stored in the ObitInfoList:
 * \li "fracPeak"OBIT_float scalar = Fraction of residual to CLEAN to 
 *                                   before major cycle [def 0.5]
 * \li "Niter"   OBIT_long scalar   = Maximum number of CLEAN iterations
 * \li "Patch"   OBIT_long scalar   = Size of beam patch in pixels [def all]
 * \li "BeamSize"OBIT_float scalar = Restoring beam FWHM (deg)
 * \li "Gain"    OBIT_float scalar = CLEAN loop gain
 * \li "minFlux" OBIT_float scalar = Minimun flux density (Jy)
 * \li "Factor"  OBIT_float array  = CLEAN depth factor
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 * \li "CCVer"   OBIT_long          = CC table version number
 * \li "noResid" OBIT_bool scalar  = If defined and TRUE, set residuals 
 *                                   to Zero in restore
 * \li "dispURL" OBIT_string scalar = URL of display server
 *
 * The OTF data file should have the following imaging parameters:
 * \li "outName" Obit_string       = Base name of image files, 
 *                                   beam = "Beam"+outname+".fits"
 *                                   clean = outname+".fits"
 * \li "outDisk" OBIT_long [1,1,1]  = Disk number for output files. [def 1]
 * \li "Scans"   OBIT_long [2,1,1]  = Range of scan numbers, [def=all]
 * \li "doCalib" OBIT_long (1,1,1)  = >0 -> calibrate,
 * \li "keepCal" OBIT_bool(1,1,1)  = True = keep cal-on data [def TRUE]
 * \li "gainUse" OBIT_long (1,1,1)  = SN/CL table version number, 0-> use highest
 * \li "flagVer" OBIT_long (1,1,1)  = Flag table version, 0-> use highest, <0-> none
 * \li "Stokes"  OBIT_string (4,1,1)= Selected output Stokes parameters:
 *                                   "I", "V", " " -> "I" [def "I"]
 * \li "BChan"   OBIT_long (1,1,1)  = First spectral channel selected. [def all]
 * \li "EChan"   OBIT_long (1,1,1)  = Highest spectral channel selected. [def all]
 * \li "Feeds"   OBIT_long (?,1,1)  = List of selected feed numbers, [def all.]
 * \li "Targets" Obit_string [?,?] = List of target names to include
 * \li "RA"      OBIT_float scalar = Center RA in deg
 * \li "Dec"     OBIT_float scalar = Center Dec in deg
 * \li "nx"      OBIT_long scalar   = Number of cells in RA
 * \li "ny"      OBIT_long scalar   = Number of cells in Declination
 * \li "beamNx"  OBIT_long scalar   = Number of "x" pixels [def 32]
 * \li "beamNy"  OBIT_long scalar   = Number of "y" pixels [def 32]
 * \li "xCells"  OBIT_float scalar = Cell spacing (deg) in RA
 * \li "yCells"  OBIT_float scalar = Cell spacing (deg) in Dec
 * \li "minWt"   OBIT_float scalar = Minimum sum of weights
 *                                   Default = 0.1.
 * \li "Proj"    OBIT_string (4,1,1) Projection string "-SIN", "-ARC", "-TAN"
 *                                   [Default "-SIN"]
 * \li "ConvType"OBIT_long scalar = Convolving function type: [def=3]
 *                  0 = pillbox, 3 = Gaussian, 4 = Exp*Sinc, 5 = Spherodial wave
 * \li "ConvParm"OBIT_float[10] = Convolving function parameters (see below)
 * 
 * Gridding convolution functions:
 * \li 0 = pillbox, 
 * \li 2 = Sinc, 
 *    Parm[0] = halfwidth in cells,
 *    Parm[1] = Expansion factor
 * \li 3 = Gaussian,
 *    Parm[0] = halfwidth in cells,[def 3.0]
 *    Parm[1] = Gaussian with as fraction or raw beam [def 1.0]
 * \li 4 = Exp*Sinc
 *    Parm[0] = halfwidth in cells, [def 2.0]
 *    Parm[1] = 1/sinc factor (cells) [def 1.55]
 *    Parm[2] = 1/exp factor (cells) [def 2.52]
 *    Parm[3] = exp power [def 2.0]
 * \li 5 = Spherodial wave
 *    Parm[0] = halfwidth in cells [def 3.0]
 *    Parm[1] = Alpha [def 5.0]
 *    Parm[2] = Expansion factor [not used]
 */

/*--------------Class definitions-------------------------------------*/
/** ObitDConCleanOTFRec Class structure. */
typedef struct {
#include "ObitDConCleanOTFRecDef.h"   /* this class definition */
} ObitDConCleanOTFRec;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitDConCleanOTF
 * returns a ObitDConCleanOTF*.
 * in = object to unreference
 */
#define ObitDConCleanOTFRecUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitDConCleanOTFRec.
 * returns a ObitDConCleanOTFRec*.
 * in = object to reference
 */
#define ObitDConCleanOTFRecRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitDConCleanOTFRecIsA(in) ObitIsA (in, ObitDConCleanOTFRecGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitDConCleanOTFRecClassInit (void);

/** Public: Default Constructor. */
ObitDConCleanOTFRec* newObitDConCleanOTFRec (gchar* name);

/** Public: Create/initialize ObitDConCleanOTFRec structures */
ObitDConCleanOTFRec* 
ObitDConCleanOTFRecCreate (gchar* name, ObitOTF *inOTF, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitDConCleanOTFRecGetClass (void);

/** Public: Copy (deep) constructor. */
ObitDConCleanOTFRec* ObitDConCleanOTFRecCopy  (ObitDConCleanOTFRec *in, 
					 ObitDConCleanOTFRec *out, ObitErr *err);

/** Public: Copy structure. */
void ObitDConCleanOTFRecClone (ObitDConCleanOTFRec *in, 
			    ObitDConCleanOTFRec *out, ObitErr *err);

/** Public: Do deconvolution. */
void ObitDConCleanOTFRecDeconvolve (ObitDCon *in, ObitErr *err);

/** Public: Get parameters. */
void  ObitDConCleanOTFRecGetParms (ObitDCon *in, ObitErr *err);

/** Public: Generate residual image */
void ObitDConCleanOTFRecSub(ObitDConClean *in, ObitErr *err);

/** Public:Select components  . */
gboolean ObitDConCleanOTFRecSelect(ObitDConClean *in, ObitErr *err);

/** Public: Subtract components and generate new residual image(s). */
void ObitDConCleanOTFRecSub(ObitDConClean *in, ObitErr *err);

/** Public:  Restore subtracted components. */
void ObitDConCleanOTFRecRestore(ObitDConClean *in, ObitErr *err);

/** Public:  Prepare for minor cycle. */
void ObitDConCleanOTFRecPixelStats(ObitDConClean *in, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitDConCleanOTFRecClassDef.h"
} ObitDConCleanOTFRecClassInfo; 

#endif /* OBITDCONCLEANOTFREC_H */ 
