/* $Id$     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2014                                          */
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
#ifndef OBITSKYMODEL_H 
#define OBITSKYMODEL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitFArray.h"
#include "ObitCArray.h"
#include "ObitCInterpolate.h"
#include "ObitImageMosaic.h"
#include "ObitUV.h"
#include "ObitTableCC.h"
#if HAVE_GPU==1  /*  GPU? */
#include "ObitGPUSkyModel.h"
#endif /* HAVE_GPU */

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSkyModel.h
 * ObitSkyModel class represents sky models and their Fourier transforms
 *
 * This class is derived from the #Obit class.
 * 
 * \section ObitSkyModelaccess Creators and Destructors
 * An ObitSkyModel will usually be created using ObitSkyModelCreate which allows 
 * specifying a name for the object as well as the ImageMosaic containing the model.
 *
 * A copy of a pointer to an ObitSkyModel should always be made using the
 * #ObitSkyModelRef function which updates the reference count in the object.
 * Then whenever freeing an ObitSkyModel or changing a pointer, the function
 * #ObitSkyModelUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 *
 * \section ObitSkyModelselect Data selection
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
 * \li  "ModelType" OBIT_int (1,1,1) Model type (ObitSkyModelType) [def OBIT_SkyModel_Comps]
 * \li  "Mode"   OBIT_int (1,1,1) Model mode (ObitSkyModelMode) [def OBIT_SkyModel_Fastest]
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
/**
 * \enum obitSkyModelType
 * enum Sky Model type, image or components, e.g. Clean components
 */
enum obitSkyModelType {
  /** Components (AIPS CC) */
  OBIT_SkyModel_Comps=0,
  /** Image */
  OBIT_SkyModel_Image,
  /** Point model */
  OBIT_SkyModel_Point
}; /* end enum obitSkyModelType */
/** typedef for enum for ObitSkyModelType. */
typedef enum obitSkyModelType ObitSkyModelType;

/**
 * \enum obitSkyModelMode
 * enum for Sky Model Component computation mode
 */
enum obitSkyModelMode {
  /** Choose fastest */
  OBIT_SkyModel_Fastest=0,
  /** DFT calculation */
  OBIT_SkyModel_DFT,  
  /** Gridded calculation */
  OBIT_SkyModel_Grid, 
  /** Some DFT and some Gridded calculation */
  OBIT_SkyModel_Mixed  
}; /* end enum obitSkyModelMode */
/** typedef for enum for ObitSkyModelMode. */
typedef enum obitSkyModelMode ObitSkyModelMode;

/**
 * \enum obitSkyModelCompType
 * enum for Sky Model Component model type
 */
enum obitSkyModelCompType {
  /** Point */
  OBIT_SkyModel_PointMod,
  /** Gaussian on sky */
  OBIT_SkyModel_GaussMod,  
  /** Convolved Gaussian */
  OBIT_SkyModel_CGaussMod,  
  /** Uniform sphere */
  OBIT_SkyModel_USphereMod, 
  /** Unknown */
  OBIT_SkyModel_Unknown,
  /** Point + spectrum */
  OBIT_SkyModel_PointModSpec = OBIT_SkyModel_PointMod+10,
  /** Gaussian on sky + spectrum */
  OBIT_SkyModel_GaussModSpec = OBIT_SkyModel_GaussMod+10,  
  /** Convolved Gaussian + spectrum */
  OBIT_SkyModel_CGaussModSpec = OBIT_SkyModel_CGaussMod+10,  
  /** Uniform sphere + spectrum */
  OBIT_SkyModel_USphereModSpec = OBIT_SkyModel_USphereMod+10 ,
  /** Point + Tabulated spectrum */
  OBIT_SkyModel_PointModTSpec = OBIT_SkyModel_PointMod+20,
  /** Gaussian on sky + Tabulated spectrum */
  OBIT_SkyModel_GaussModTSpec = OBIT_SkyModel_GaussMod+20,  
  /** Convolved Gaussian + sTabulated pectrum */
  OBIT_SkyModel_CGaussModTSpec = OBIT_SkyModel_CGaussMod+20,  
  /** Uniform sphere + Tabulated spectrum */
  OBIT_SkyModel_USphereModTSpec = OBIT_SkyModel_USphereMod+20 
}; /* end enum obitSkyModelCompType */
/** typedef for enum for ObitSkyModelCompType. */
typedef enum obitSkyModelCompType ObitSkyModelCompType;

/*--------------Class definitions-------------------------------------*/
/** ObitSkyModel Class structure. */
typedef struct {
#include "ObitSkyModelDef.h"   /* this class definition */
} ObitSkyModel;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitSkyModel
 * returns a ObitSkyModel*.
 * in = object to unreference
 */
#define ObitSkyModelUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitSkyModel.
 * returns a ObitSkyModel*.
 * in = object to reference
 */
#define ObitSkyModelRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitSkyModelIsA(in) ObitIsA (in, ObitSkyModelGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitSkyModelClassInit (void);

/** Public: Default Constructor. */
ObitSkyModel* newObitSkyModel (gchar* name);

/** Public: Create SkyModel object from description in an ObitInfoList */
ObitSkyModel* ObitSkyModelFromInfo (gchar *prefix, ObitInfoList *inList, 
				    ObitErr *err);
typedef ObitSkyModel* (*ObitSkyModelFromInfoFP) (gchar *prefix, 
						 ObitInfoList *inList, 
						 ObitErr *err);

/** Public: Create/initialize ObitSkyModel structures */
ObitSkyModel* ObitSkyModelCreate (gchar* name, ObitImageMosaic* mosaic);
/** Typedef for definition of class pointer structure */
typedef ObitSkyModel* (*ObitSkyModelCreateFP) (gchar* name, ObitImageMosaic* mosaic);

/** Public: initialize ObitSkyModel structures */
void ObitSkyModelInitMod (ObitSkyModel* in, ObitUV *uvdata, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitSkyModelInitModFP) (ObitSkyModel* in, ObitUV *uvdata, 
				       ObitErr *err);

/** Public: shutdown ObitSkyModel processes */
void ObitSkyModelShutDownMod (ObitSkyModel* in, ObitUV *uvdata, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitSkyModelShutDownModFP) (ObitSkyModel* in, ObitUV *uvdata, 
					   ObitErr *err);

/** Public: initialize model for pass in time through data */
void ObitSkyModelInitModel (ObitSkyModel* in, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitSkyModelInitModelFP) (ObitSkyModel* in, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitSkyModelGetClass (void);

/** Public: Copy (deep) constructor. */
ObitSkyModel* ObitSkyModelCopy  (ObitSkyModel *in, ObitSkyModel *out, ObitErr *err);

/** Public: Copy structure. */
void ObitSkyModelClone (ObitSkyModel *in, ObitSkyModel *out, ObitErr *err);

/** Public: Subtract model from an ObitUV */
ObitIOCode ObitSkyModelSubUV (ObitSkyModel *in, ObitUV *indata, ObitUV *outdata, 
			      ObitErr *err);
typedef ObitIOCode (*ObitSkyModelSubUVFP) (ObitSkyModel *in, ObitUV *indata, 
					  ObitUV *outdata, ObitErr *err);

/** Public: Divide model into an ObitUV */
ObitIOCode ObitSkyModelDivUV (ObitSkyModel *in, ObitUV *indata, ObitUV *outdata, 
			      ObitErr *err);
typedef ObitIOCode (*ObitSkyModelDivUVFP) (ObitSkyModel *in, ObitUV *indata, 
					  ObitUV *outdata, ObitErr *err);

/** Public: Load specified image and plane */
gboolean ObitSkyModelLoad (ObitSkyModel *in, olong image, ObitUV *uvdata,
			     ObitErr *err);
typedef gboolean (*ObitSkyModelLoadFP) (ObitSkyModel *in, olong field,
					ObitUV *uvdata, ObitErr *err);

/** Public: Calculate Fourier transform of model for current buffer in an ObitUV */
void ObitSkyModelFT (ObitSkyModel *in, olong field,ObitUV *uvdata, ObitErr *err);
typedef void (*ObitSkyModelFTFP) (ObitSkyModel *in, olong field,ObitUV *uvdata, 
				  ObitErr *err);

/** Public: Load point model, may be overridden in derived class */
gboolean ObitSkyModelLoadPoint (ObitSkyModel *in, ObitUV *uvdata, ObitErr *err);
typedef gboolean (*ObitSkyModelLoadPointFP) (ObitSkyModel *in, ObitUV *uvdata, 
					     ObitErr *err);

/** Public: Load  Components model, may be overridden in derived class */
gboolean ObitSkyModelLoadComps (ObitSkyModel *in, olong n, ObitUV *uvdata, 
				  ObitErr *err);
typedef gboolean (*ObitSkyModelLoadCompsFP) (ObitSkyModel *in, olong n, 
					      ObitUV *uvdata, ObitErr *err);

/** Public: Grid  Components model, may be overridden in derived class */
gboolean ObitSkyModelGridComps (ObitSkyModel *in, olong n, ObitUV *uvdata, 
				  ObitErr *err);
typedef gboolean (*ObitSkyModelGridCompsFP) (ObitSkyModel *in, olong n, 
					       ObitUV *uvdata, ObitErr *err);

/** Public: Load/FT image model, may be overridden in derived class */
gboolean ObitSkyModelLoadImage (ObitSkyModel *in, olong n, ObitUV *uvdata, 
				  ObitErr *err);
typedef gboolean (*ObitSkyModelLoadImageFP) (ObitSkyModel *in, olong n, 
					      ObitUV *uvdata, ObitErr *err);

/** Public: FT by DFT, may be overridden in derived class */
void ObitSkyModelFTDFT (ObitSkyModel *in, olong field, ObitUV *uvdata, ObitErr *err);
typedef void (*ObitSkyModelFTDFTFP) (ObitSkyModel *in, olong field, ObitUV *uvdata, 
					  ObitErr *err);

/** Public: FT by Gridding, may be overridden in derived class */
void ObitSkyModelFTGrid (ObitSkyModel *in, olong field, ObitUV *uvdata, ObitErr *err);
typedef void (*ObitSkyModelFTGridFP) (ObitSkyModel *in, olong field, ObitUV *uvdata, 
					   ObitErr *err);

/** Public: Sum flux in Clean Model */
ofloat ObitSkyModelSum (ObitSkyModel *in, ObitErr *err);
typedef ofloat (*ObitSkyModelSumFP) (ObitSkyModel *in, ObitErr *err);

/** Public: Compress CC Tables */
void ObitSkyModelCompressCC (ObitSkyModel *in, ObitErr *err);
typedef void (*ObitSkyModelCompressCCFP) (ObitSkyModel *in, ObitErr *err);

/** Public: Get input parameters from info */
void  ObitSkyModelGetInput (ObitSkyModel* in, ObitErr *err);
typedef void (*ObitSkyModelGetInputFP) (ObitSkyModel* in, ObitErr *err);

/** Public: Decide model method */
void  ObitSkyModelChose (ObitSkyModel* in, ObitUV* uvdata);
typedef void (*ObitSkyModelChoseFP) (ObitSkyModel* in, ObitUV* uvdata);

/** Public: Fill in data selection values */
void ObitSkyModelSetSelect (ObitSkyModel* in, ObitUV* uvdata, ObitErr *err);
typedef void (*ObitSkyModelSetSelectFP) (ObitSkyModel* in, ObitUV* uvdata, ObitErr *err);

/** Public: Decide next block of channels if doing PB correction */
gboolean ObitSkyModelsetPBChans(ObitSkyModel* in, ObitUV* uvdata, ObitErr *err);
typedef gboolean (*ObitSkyModelsetPBChansFP) (ObitSkyModel* in, ObitUV* uvdata, ObitErr *err);

/** Public: return ObitTableCC with possible PB corrections */
ObitTableCC* ObitSkyModelgetPBCCTab (ObitSkyModel* in, ObitUV* uvdata, 
				     olong field, olong *inCCVer, olong *outCCver,
				     olong *startCC, olong *endCC, ofloat range[2],
				     ObitErr *err);
typedef ObitTableCC* 
(*ObitSkyModelgetPBCCTabFP) (ObitSkyModel* in, ObitUV* uvdata, 
			     olong field, olong *inCCVer, olong *outCCver,
			     olong *startCC, olong *endCC, ObitErr *err);

/** Public: fill in->plane with image and possibly PB corrected */
void ObitSkyModelgetPBImage (ObitSkyModel* in, ObitUV* uvdata, olong field, 
			     ObitErr *err);
typedef void (*ObitSkyModelgetPBImageFP) (ObitSkyModel* in, ObitUV* uvdata, olong field, 
					  ObitErr *err);

/** Public: Grid/FT components */
gboolean ObitSkyModelGridFTComps (ObitSkyModel* in, olong field, ObitUV* uvdata, 
				  ObitErr *err);
typedef gboolean (*ObitSkyModelGridFTCompsFP) (ObitSkyModel* in, olong field, ObitUV* uvdata, 
				  ObitErr *err);

/** Public: Load Grid components */
void ObitSkyModelLoadGridComps (ObitSkyModel* in, olong field, ObitUV* uvdata, 
				ofloat gparm[3], olong *ncomp, ObitErr *err);
typedef void (*ObitSkyModelLoadGridCompsFP) (ObitSkyModel* in, olong field, ObitUV* uvdata, 
					     ofloat gparm[3], olong *ncomp, ObitErr *err);

/** Public: FT image array in in->plane */
void ObitSkyModelFTImage (ObitSkyModel* in, ObitFArray *inArray, 
			  ObitCArray *outArray);
typedef void (*ObitSkyModelFTImageFP) (ObitSkyModel* in, ObitFArray *inArray, 
				       ObitCArray *outArray);

/** Public: Add field to mosaic */
void ObitSkyModelAddField (ObitSkyModel *in, ObitErr *err);
typedef void (*ObitSkyModelAddFieldFP) (ObitSkyModel* in, ObitErr *err);

/** Public: Extract information about underlying structures to ObitInfoList */
void ObitSkyModelGetInfo (ObitSkyModel *in, gchar *prefix, ObitInfoList *outList, 
			  ObitErr *err);
typedef void 
(*ObitSkyModelGetInfoFP) (ObitSkyModel *in, gchar *prefix, ObitInfoList *outList, 
			  ObitErr *err);
/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitSkyModelClassDef.h"
} ObitSkyModelClassInfo; 

#endif /* OBITFSKYMODEL_H */ 
