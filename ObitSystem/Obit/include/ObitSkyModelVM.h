/* $Id$     */
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
#ifndef OBITSKYMODELVM_H 
#define OBITSKYMODELVM_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitFArray.h"
#include "ObitCArray.h"
#include "ObitCInterpolate.h"
#include "ObitImageMosaic.h"
#include "ObitUV.h"
#include "ObitSkyModel.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSkyModelVM.h
 * ObitSkyModelVMIon class represents temporally or spatially variable effect 
 * corrected sky models and their Fourier transforms .
 *
 * This is a virtual ObitSkyModel class for implementing incorporation of 
 * temporally and/or spatially variable effects.
 *
 * This class is derived from the #ObitSkyModel class.
 *
 * \section ObitSkyModelVMaccess Creators and Destructors
 * An ObitSkyModelVM will usually be created using ObitSkyModelVMCreate which allows 
 * specifying a name for the object as well as the ImageMosaic containing the model.
 *
 * A copy of a pointer to an ObitSkyModelVM should always be made using the
 * #ObitSkyModelVMRef function which updates the reference count in the object.
 * Then whenever freeing an ObitSkyModelVM or changing a pointer, the function
 * #ObitSkyModelVMUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 *
 * \section ObitSkyModelVMselect Data selection
 * The selection of data to be modified is through values added to the info 
 *(#ObitInfoList) member and include the following:
 * \li  "Stokes" OBIT_string (4,1,1) Selected output Stokes parameters:
 *               "    "=> no translation,"I   ","V   ","Q   ", "U   ", 
 *               "IQU ", "IQUV",  "IV  ", "RR  ", "LL  ", "RL  ", "LR  ", 
 *               "HALF" = RR,LL, "FULL"=RR,LL,RL,LR. [default "I"]
 *               In the above 'F' can substitute for "formal" 'I' (both RR+LL).
 * \li  "ionVer"OBIT_int (1,1,1) NI table version number, 0-> use highest, def=1
 * \li  "BChan" OBIT_int (1,1,1) First spectral channel selected. [def all]
 * \li  "EChan" OBIT_int (1,1,1) Highest spectral channel selected. [def all]
 * \li  "BIF"   OBIT_int (1,1,1) First "IF" selected. [def all]
 * \li  "EIF"   OBIT_int (1,1,1) Highest "IF" selected. [def all]
 * \li  "do3D"  OBIT_bool (1,1,1) If 3D imaging wanted. [def false]
 * \li  "DIVIDE" OBIT_bool (1,1,1) If division rather than subtraction wanted. [def false]
 * \li  "REPLACE" OBIT_bool (1,1,1) If TRUE replace data with model. [def false]
 * \li  "PBCor"  OBIT_bool (1,1,1) If TRUE make relative primary beam corrections [def true]
 * \li  "noNeg"  OBIT_bool (1,1,1) If TRUE only positive flux comps are to be used [def false]
 * \li  "antSize" OBIT_float (1,1,1) Diameter of antennas for rel. r,.[def 25.0]
 * \li  "Factor" OBIT_float (1,1,1) model multiplications factor (-1=>add) [def 1]
 * \li  "minFlux" OBIT_float (1,1,1) Minimum flux density model or pixel [def -1.0e20]
 * \li  "ModelType" OBIT_int (1,1,1) Model type (ObitSkyModelType) [def OBIT_SkyModel_Comps]
 * \li  "MODPTFLX" OBIT_float (1,1,1) Point model flux in Jy, [def 0.0]
 * \li  "MODPTXOF" OBIT_float (1,1,1) Point model "x" offset in deg  [def 0.0]
 * \li  "MODPTYOF" OBIT_float (1,1,1) Point model "y" offset in deg  [def 0.0]
 * \li  "MODPTYPM" OBIT_float (4,1,1) Point other parameters  [def all 0.0]
 * \li  "CCVer" OBIT_int (?,1,1) CC table versions to use [def all 0 => highest]
 * \li  "BComp" OBIT_int (?,1,1) Start CC to use per table, 1-rel [def 1 ]
 * \li  "EComp" OBIT_int (?,1,1) Highest CC to use per table, 1-rel [def to end ]
 * \li  "UpdateInt" OBIT_float (1,1,1) Model update interval (min)  [def 1 min]
 */

/*-------------- enumerations -------------------------------------*/

/*--------------Class definitions-------------------------------------*/
/** ObitSkyModelVM Class structure. */
typedef struct {
#include "ObitSkyModelVMDef.h"   /* this class definition */
} ObitSkyModelVM;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitSkyModelVM
 * returns a ObitSkyModelVM*.
 * in = object to unreference
 */
#define ObitSkyModelVMUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitSkyModelVM.
 * returns a ObitSkyModelVM*.
 * in = object to reference
 */
#define ObitSkyModelVMRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitSkyModelVMIsA(in) ObitIsA (in, ObitSkyModelVMGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitSkyModelVMClassInit (void);

/** Public: Default Constructor. */
ObitSkyModelVM* newObitSkyModelVM (gchar* name);

/** Public: Create/initialize ObitSkyModelVM structures */
ObitSkyModelVM* ObitSkyModelVMCreate (gchar* name, ObitImageMosaic* mosaic);

/** Public: Get Inputs. */
void  ObitSkyModelVMGetInput (ObitSkyModel* inn, ObitErr *err);

/** Public: initialize ObitSkyModelVM structures */
void ObitSkyModelVMInitMod (ObitSkyModel* in, ObitUV *uvdata, ObitErr *err);

/** Public: shutdown ObitSkyModel processes */
void ObitSkyModelVMShutDownMod (ObitSkyModel* in,ObitUV *uvdata, ObitErr *err);

/** Public: initialize model for pass in time through data */
void ObitSkyModelVMInitModel (ObitSkyModel* in, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitSkyModelVMGetClass (void);

/** Public: Copy (deep) constructor. */
ObitSkyModel* ObitSkyModelVMCopy  (ObitSkyModel *in, ObitSkyModel *out, ObitErr *err);

/** Public: Routine to update model */
void ObitSkyModelVMUpdateModel (ObitSkyModelVM *in, ofloat time, olong suba,
				ObitUV *uvdata, ObitErr *err);
typedef void (*ObitSkyModelVMUpdateModelFP) (ObitSkyModelVM *in, 
					     ofloat time, olong suba,
					     ObitUV *uvdata, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitSkyModelVMClassDef.h"
} ObitSkyModelVMClassInfo; 

#endif /* OBITFSKYMODELVM_H */ 
