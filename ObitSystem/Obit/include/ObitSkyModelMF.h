/* $Id: ObitSkyModelMF.h 61 2008-12-19 18:14:49Z bill.cotton $     */
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
#ifndef OBITSKYMODELMF_H 
#define OBITSKYMODELMF_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitFArray.h"
#include "ObitCArray.h"
#include "ObitCInterpolate.h"
#include "ObitImageMosaic.h"
#include "ObitUV.h"
#include "ObitTableCC.h"
#include "ObitSkyModel.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSkyModelMF.h
 * ObitSkyModelMF class represents sky models and their Fourier transforms
 *
 * This class is derived from the #Obit class.
 * 
 * \section ObitSkyModelMFaccess Creators and Destructors
 * An ObitSkyModelMF will usually be created using ObitSkyModelMFCreate which allows 
 * specifying a name for the object as well as the ImageMosaic containing the model.
 *
 * A copy of a pointer to an ObitSkyModelMF should always be made using the
 * #ObitSkyModelMFRef function which updates the reference count in the object.
 * Then whenever freeing an ObitSkyModelMF or changing a pointer, the function
 * #ObitSkyModelMFUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 *
 * \section ObitSkyModelMFselect Data selection
 * The selection of data to be modified is through values added to the info 
 *(#ObitInfoList) member and include the following:
 * \li  "Stokes" OBIT_string (4,1,1) Selected output Stokes parameters:
 *               "    "=> no translation,"I   ","V   ","Q   ", "U   ", 
 *               "IQU ", "IQUV",  "IV  ", "RR  ", "LL  ", "RL  ", "LR  ", 
 *               "HALF" = RR,LL, "FULL"=RR,LL,RL,LR. [default "I"]
 *               In the above 'F' can substitute for "formal" 'I' (both RR+LL).
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
 * \li  "ModelType" OBIT_int (1,1,1) Model type (ObitSkyModelMFType) [def OBIT_SkyModelMF_Comps]
 * \li  "Mode"   OBIT_int (1,1,1) Model mode (ObitSkyModelMFMode) [def OBIT_SkyModelMF_Fastest]
 * \li  "MODPTFLX" OBIT_float (1,1,1) Point model flux in Jy, [def 0.0]
 * \li  "MODPTXOF" OBIT_float (1,1,1) Point model "x" offset in deg  [def 0.0]
 * \li  "MODPTYOF" OBIT_float (1,1,1) Point model "y" offset in deg  [def 0.0]
 * \li  "MODPTYPM" OBIT_float (4,1,1) Point other parameters  [def all 0.0]
 * \li  "CCVer" OBIT_int (?,1,1) CC table versions to use [def all 0 => highest]
 * \li  "BComp" OBIT_int (?,1,1) Start CC to use per table, 1-rel [def 1 ]
 * \li  "EComp" OBIT_int (?,1,1) Highest CC to use per table, 1-rel [def to end ]
 * \li "prtLv" OBIT_int message level  [def 0]
 */

/*-------------- enumerations -------------------------------------*/

/*--------------Class definitions-------------------------------------*/
/** ObitSkyModelMF Class structure. */
typedef struct {
#include "ObitSkyModelMFDef.h"   /* this class definition */
} ObitSkyModelMF;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitSkyModelMF
 * returns a ObitSkyModelMF*.
 * in = object to unreference
 */
#define ObitSkyModelMFUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitSkyModelMF.
 * returns a ObitSkyModelMF*.
 * in = object to reference
 */
#define ObitSkyModelMFRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitSkyModelMFIsA(in) ObitIsA (in, ObitSkyModelMFGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitSkyModelMFClassInit (void);

/** Public: Default Constructor. */
ObitSkyModelMF* newObitSkyModelMF (gchar* name);

/** Public: Create SkyModelMF object from description in an ObitInfoList */
ObitSkyModelMF* ObitSkyModelMFFromInfo (gchar *prefix, ObitInfoList *inList, 
					ObitErr *err);

/** Public: Create/initialize ObitSkyModelMF structures */
ObitSkyModelMF* ObitSkyModelMFCreate (gchar* name, ObitImageMosaic* mosaic);

/** Public: initialize ObitSkyModelMF structures */
void ObitSkyModelMFInitMod (ObitSkyModel* in, ObitUV *uvdata, ObitErr *err);

/** Public: shutdown ObitSkyModelMF processes */
void ObitSkyModelMFShutDownMod (ObitSkyModel* in, ObitUV *uvdata, ObitErr *err);

/** Public: initialize model for pass in time through data */
void ObitSkyModelMFInitModel (ObitSkyModel* in, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitSkyModelMFGetClass (void);

/** Public: Copy (deep) constructor. */
ObitSkyModel* ObitSkyModelMFCopy  (ObitSkyModel *in, ObitSkyModel *out, ObitErr *err);

/** Public: Copy structure. */
void ObitSkyModelMFClone (ObitSkyModel *in, ObitSkyModel *out, ObitErr *err);

/** Public: Subtract model from an ObitUV */
ObitIOCode ObitSkyModelMFSubUV (ObitSkyModel *in, ObitUV *indata, ObitUV *outdata, 
			      ObitErr *err);

/** Public: Divide model into an ObitUV */
ObitIOCode ObitSkyModelMFDivUV (ObitSkyModel *in, ObitUV *indata, ObitUV *outdata, 
			      ObitErr *err);

/** Public: Load specified image and plane */
gboolean ObitSkyModelMFLoad (ObitSkyModel *in, olong image, ObitUV *uvdata,
			     ObitErr *err);

/** Public: Load point model, may be overridden in derived class */
gboolean ObitSkyModelMFLoadPoint (ObitSkyModel *in, ObitUV *uvdata, ObitErr *err);

/** Public: Load  Components model, may be overridden in derived class */
gboolean ObitSkyModelMFLoadComps (ObitSkyModel *in, olong n, ObitUV *uvdata, 
				  ObitErr *err);

/** Public: Grid  Components model, may be overridden in derived class */
gboolean ObitSkyModelMFGridComps (ObitSkyModel *in, olong n, ObitUV *uvdata, 
				  ObitErr *err);

/** Public: Load/FT image model, may be overridden in derived class */
gboolean ObitSkyModelMFLoadImage (ObitSkyModel *in, olong n, ObitUV *uvdata, 
				  ObitErr *err);

/** Public: FT by DFT, may be overridden in derived class */
void ObitSkyModelMFFTDFT (ObitSkyModel *in, olong field, ObitUV *uvdata, ObitErr *err);

/** Public: FT by Gridding, may be overridden in derived class */
void ObitSkyModelMFFTGrid (ObitSkyModel *in, olong field, ObitUV *uvdata, ObitErr *err);

/** Public: Sum flux in Clean Model */
ofloat ObitSkyModelMFSum (ObitSkyModel *in, ObitErr *err);

/** Public: Get input parameters from info */
void  ObitSkyModelMFGetInput (ObitSkyModel* in, ObitErr *err);

/** Public: Decide model method */
void  ObitSkyModelMFChose (ObitSkyModel* in, ObitUV* uvdata);

/** Public: Fill in data selection values */
void ObitSkyModelMFSetSelect (ObitSkyModel* in, ObitUV* uvdata, ObitErr *err);

/** Public: fill in->plane with image and possibly PB corrected */
void ObitSkyModelMFgetPBImage (ObitSkyModel* in, ObitUV* uvdata, olong field, 
			     ObitErr *err);

/** Public: Grid/FT components */
gboolean ObitSkyModelMFGridFTComps (ObitSkyModel* in, olong field, ObitUV* uvdata, 
				  ObitErr *err);

/** Public: Load Grid components */
void ObitSkyModelMFLoadGridComps (ObitSkyModel* in, olong field, ObitUV* uvdata, 
				ofloat gparm[3], olong *ncomp, ObitErr *err);

/** Public: FT image array in in->plane */
void ObitSkyModelMFFTImage (ObitSkyModel* in, ObitFArray *inArray, 
			  ObitCArray *outArray);

/** Public: Add field to mosaic */
void ObitSkyModelMFAddField (ObitSkyModel *in, ObitErr *err);

/** Public: Extract information about underlying structures to ObitInfoList */
void ObitSkyModelMFGetInfo (ObitSkyModel *in, gchar *prefix, ObitInfoList *outList, 
			  ObitErr *err);
/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitSkyModelMFClassDef.h"
} ObitSkyModelMFClassInfo; 

#endif /* OBITFSKYMODELMF_H */ 
