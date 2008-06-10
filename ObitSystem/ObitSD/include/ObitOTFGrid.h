/* $Id: ObitOTFGrid.h,v 1.14 2008/05/15 21:31:49 bcotton Exp $  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
#ifndef OBITOTFGRID_H 
#define OBITOTFGRID_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitOTF.h"
#include "ObitFArray.h"
#include "ObitImage.h"
#include "ObitOTFSkyModel.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFGrid.h
 * ObitOTFGrid uv data class definition.
 *
 * This class is for creating and manipulating a OTFGrid as a memory resident 
 * intermediate entity between ObitOTF and ObitImage classes.
 * The beam corresponding to each image should be made first using the
 * same ObitOTFGrid.
 * 
 * \section ObitOTFGridaccess Creators and Destructors
 * An ObitOTFGrid can be created using newObitOTFGrid which allows specifying 
 * a name for the object.  This name is used to label messages.
 *
 * A copy of a pointer to an ObitOTFGrid should always be made using the
 * #ObitOTFGridRef function which updates the reference count in the object.
 * Then whenever freeing an ObitOTFGrid or changing a pointer, the function
 * #ObitOTFGridUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 *
 * \section ObitOTFGridparameters Control Parameters
 * The imaging control parameters are passed through the info object 
 * on the OTF data, these control both the output image files and the 
 * processing parameters.
 * OTF Data selection/calibration/editing control
 * \li  "doCalSelect" OBIT_bool (1,1,1) Select/calibrate/edit data?
 * \li  "doCalib" OBIT_int (1,1,1) >0 -> calibrate,
 * \li  "gainUse" OBIT_int (1,1,1) SN/CL table version number, 0-> use highest
 * \li  "flagVer" OBIT_int (1,1,1) Flag table version, 0-> use highest, <0-> none
 * \li  "Stokes" OBIT_string (4,1,1) Selected output Stokes parameters:
 *               "I", "V", " " -> "I"
 * \li  "BChan" OBIT_int (1,1,1) First spectral channel selected. [def all]
 * \li  "EChan" OBIT_int (1,1,1) Highest spectral channel selected. [def all]
 * \li  "Targets" OBIT_string (?,?,1) Target names selected. [def all]
 * \li  "timeRange" OBIT_float (2,1,1) Selected timerange in days. [def all]
 * \li  "Scans" OBIT_int (2,1,1) Lowest and highest selected scan numbers. [def all]
 * \li  "Feeds" OBIT_int (?,1,1) a list of selected feed numbers, [def all.]
 * 
 * Gridding/imaging parameters
 * \li "ConvType"   OBIT_long scalar = Convolving function type: [def=3]
 *                  0 = pillbox, 3 = Gaussian, 4 = Exp*Sinc, 5 = Spherodial wave
 * \li "ConvParm" OBIT_float [10] = Convolving function parameters
 *                  of the antenna beam size, Default = 1.0.
 * \li "minWt"     OBIT_float (1,1,1) Minimum summed gridding convolution weight 
 *                 as a fraction of the maximum [def 0.01]
 * \li "Clip"      OBIT_float scalar = data values with abs. value larger are set zero weight
 * \li "beamNx"    OBIT_int scalar Number of "x" pixels [def 32]
 * \li "beamNy"    OBIT_int scalar Number of "y" pixels [def 32]
 * \li "doScale"   OBIT_bool scalar If true, convolve/scale beam [def TRUE]
 *
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
/** ObitOTFGrid Class structure. */
typedef struct {
#include "ObitOTFGridDef.h"   /* this class definition */
} ObitOTFGrid;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitOTFGrid
 * returns a ObitOTFGrid*.
 * in = object to unreference
 */
#define ObitOTFGridUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitOTFGrid.
 * returns a ObitOTFGrid*.
 * in = object to reference
 */
#define ObitOTFGridRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitOTFGridIsA(in) ObitIsA (in, ObitOTFGridGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitOTFGridClassInit (void);

/** Public: Constructor. */
ObitOTFGrid* newObitOTFGrid (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitOTFGridGetClass (void);

/** Public: initialize/reset ObitOTFGrid structures */
void ObitOTFGridSetup (ObitOTFGrid *in, ObitOTF *OTFin, 
		       ObitImageDesc *imageDesc, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitOTFGridSetupFP) (ObitOTFGrid *in, ObitOTF *OTFin, 
				    ObitImageDesc *imageDesc, ObitErr *err);

/** Public: Read uv data accumulating to grid */
void ObitOTFGridReadOTF (ObitOTFGrid *in, ObitOTF *OTFin, ObitErr *err);
typedef ObitIOCode (*ObitOTFGridReadOTFFP) (ObitOTFGrid *in, ObitOTF *OTFin, 
					    ObitErr *err);

/** Public: Normalize image grid with weight grid. */
void ObitOTFGridNorm (ObitOTFGrid *in, ObitFArray *array, ObitImageDesc *imageDesc, 
		      ObitErr *err);
typedef void (*ObitOTFGridNormFP) (ObitOTFGrid *in, ObitFArray *array, 
				   ObitImageDesc *imageDesc, ObitErr *err);

/** Public: Make PSF ("Beam"). */
void ObitOTFGridMakeBeam (ObitOTFGrid* in, ObitImage *image, ObitImage *Beam, 
			  ObitErr *err);
typedef void (*ObitOTFGridMakeBeamFP) (ObitOTFGrid* in, ObitImage *image, 
				       ObitImage *Beam, ObitErr *err);
/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitOTFGridClassDef.h"
} ObitOTFGridClassInfo; 

#endif /* OBITOTFGRID_H */ 
