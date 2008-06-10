/* $Id$  */
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
#ifndef OBITSKYMODELVMION_H 
#define OBITSKYMODELVMION_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitFArray.h"
#include "ObitCArray.h"
#include "ObitCInterpolate.h"
#include "ObitImageMosaic.h"
#include "ObitUV.h"
#include "ObitSkyModelVM.h"
#include "ObitTableNI.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSkyModelVMIon.h
 * ObitSkyModelVMIon class represents ionospheric corrected sky models and their Fourier transforms 
 *
 * This class represents sky models incorporating ionospheric phase models and 
 * their Fourier transforms and is derived from the ObitSkyModelVMIon class.
 * Models of visibility data incorporating the corrupting effects of a 
 * time and spatially variable ionospheric model are calculated
 *
 * This class is derived from the #ObitSkyModelVM class.
 *
 * \section ObitSkyModelVMIonaccess Creators and Destructors
 * An ObitSkyModelVMIon will usually be created using ObitSkyModelVMIonCreate which allows 
 * specifying a name for the object as well as the ImageMosaic containing the model.
 *
 * A copy of a pointer to an ObitSkyModelVMIon should always be made using the
 * #ObitSkyModelVMIonRef function which updates the reference count in the object.
 * Then whenever freeing an ObitSkyModelVMIon or changing a pointer, the function
 * #ObitSkyModelVMIonUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 *
 * \section ObitSkyModelVMIonselect Data selection
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
 * \li  "ionVer"OBIT_int (1,1,1) NI table version number, 0-> use highest, def=1
 */

/*-------------- enumerations -------------------------------------*/

/*--------------Class definitions-------------------------------------*/
/** ObitSkyModelVMIon Class structure. */
typedef struct {
#include "ObitSkyModelVMIonDef.h"   /* this class definition */
} ObitSkyModelVMIon;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitSkyModelVMIon
 * returns a ObitSkyModelVMIon*.
 * in = object to unreference
 */
#define ObitSkyModelVMIonUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitSkyModelVMIon.
 * returns a ObitSkyModelVMIon*.
 * in = object to reference
 */
#define ObitSkyModelVMIonRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitSkyModelVMIonIsA(in) ObitIsA (in, ObitSkyModelVMIonGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitSkyModelVMIonClassInit (void);

/** Public: Default Constructor. */
ObitSkyModelVMIon* newObitSkyModelVMIon (gchar* name);

/** Public: Create/initialize ObitSkyModelVMIon structures */
ObitSkyModelVMIon* ObitSkyModelVMIonCreate (gchar* name, ObitImageMosaic* mosaic);

/** Public: initialize ObitSkyModelVMIon structures */
void ObitSkyModelVMIonInitMod (ObitSkyModel* in, ObitUV *uvdata, ObitErr *err);

/** Public: shutdown ObitSkyModel processes */
void ObitSkyModelVMIonShutDownMod (ObitSkyModel* in, ObitUV *uvdata, ObitErr *err);

/** Public: initialize model for pass through data */
void ObitSkyModelVMInitModel (ObitSkyModel* in, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitSkyModelVMIonGetClass (void);

/** Public: Copy (deep) constructor. */
ObitSkyModelVMIon* ObitSkyModelVMIonCopy  (ObitSkyModelVMIon *in, ObitSkyModelVMIon *out, ObitErr *err);

/** Public: Copy structure. */
void ObitSkyModelVMIonClone (ObitSkyModelVMIon *in, ObitSkyModelVMIon *out, 
			     ObitErr *err);

/** Public: Routine to update model */
void ObitSkyModelVMIonUpdateModel (ObitSkyModelVM *in, ofloat time, olong suba,
				   ObitUV *uvdata, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitSkyModelVMIonClassDef.h"
} ObitSkyModelVMIonClassInfo; 

#endif /* OBITFSKYMODELVMION_H */ 
