/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009-2024                                          */
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
#ifndef OBITSKYMODELVMBEAM_H 
#define OBITSKYMODELVMBEAM_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitFArray.h"
#include "ObitCArray.h"
#include "ObitCInterpolate.h"
#include "ObitImageMosaic.h"
#include "ObitUV.h"
#include "ObitSkyModelVM.h"
#include "ObitImageInterp.h"
#include "ObitBeamShape.h"
#include "ObitMatx.h"
#include "ObitComplex.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSkyModelVMBeam.h
 * ObitSkyModelVMIon class represents VLA beam squint corrected sky models 
 * and their Fourier transforms 
 *
 * This class represents sky models incorporating full beam polarization corrections
 * for NRAO-like antennas with offset feeds.  
 *
 * This class is derived from the #ObitSkyModelVM class.
 *
 * The accurate model calculation uses a Direct Fourier Transform (DFT) to
 * calculate the effects of the beam gain and polarization effects based on
 * images of these effects as a function of frequency.
 * As such, this class can only operate on data sets containing RR, LL, RL and LR
 * correlations.  Requested Stokes should always be "    ".
 * Since these effects limit dynamic range, the correction can be applied 
 * only to the brighter parts of a given CLEAN model and faster methods may
 * be applied to the less bright parts of the model.
 * Control over when the accurate, slow and less accurate but faster methods 
 * are used is controlled by InfoList members Threshold and maxResid (see below).
 * Threshold sets the level for the accurate model.    There are several cases:
 * \li CLEANing
 * For use in CLEANing, maxResid should be either unspecified or set to a negative 
 * number.  In this case, the value used will be derived from the residual images.
 * Whenever the maximum residual level at the start of a CLEAN major cycle exceeds
 * Threshold, all component subtraction will use the high accuracy model.
 * Otherwise, the modeling type is controlled by the info item Mode.
 * \li Calculation of complete model.
 * If a complete model calculation such as for a self calibration or model subtraction
 * is needed then explicitly set maxResid to 0.0.  In this case, the component lists
 * will be compressed (all components in a given cell combined) and the decision to use
 * the accurate or fast model calculation is made based on whether the component has 
 * an absolute value of flux density above or below Threshold.
 *
 * \section ObitSkyModelVMBeamaccess Creators and Destructors
 * An ObitSkyModelVMBeam will usually be created using ObitSkyModelVMBeamCreate 
 * which allows specifying a name for the object as well as the ImageMosaic 
 * containing the model.
 *
 * A copy of a pointer to an ObitSkyModelVMBeam should always be made using the
 * #ObitSkyModelVMBeamRef function which updates the reference count in the object.
 * Then whenever freeing an ObitSkyModelVMBeam or changing a pointer, the function
 * #ObitSkyModelVMBeamUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 *
 * \section ObitSkyModelVMBeamselect Data selection
 * The selection of data to be modified is through values added to the info 
 *(#ObitInfoList) member and include the following:
 * \li  "Stokes" OBIT_string (4,1,1) Selected output Stokes parameters:
 *               "    "=> no translation, [default "    "]
 * \li  "ionVer"OBIT_int (1,1,1) NI table version number, 0-> use highest, def=1
 * \li  "BChan" OBIT_int (1,1,1) First spectral channel selected. [def all]
 * \li  "EChan" OBIT_int (1,1,1) Highest spectral channel selected. [def all]
 * \li  "BIF"   OBIT_int (1,1,1) First "IF" selected. [def all]
 * \li  "EIF"   OBIT_int (1,1,1) Highest "IF" selected. [def all]
 * \li  "do3D"  OBIT_bool (1,1,1) If 3D imaging wanted. [def false]
 * \li  "DIVIDE" OBIT_bool (1,1,1) If division rather than subtraction wanted. [def false]
 * \li  "REPLACE" OBIT_bool (1,1,1) If TRUE replace data with model. [def false]
 * \li  "PBCor"   OBIT_bool (1,1,1) If TRUE make relative primary beam corrections [def true]
 * \li  "noNeg"   OBIT_bool (1,1,1) If TRUE only positive flux comps are to be used [def false]
 * \li  "antSize" OBIT_float (1,1,1) Diameter of antennas for rel. r,.[def 25.0]
 * \li  "PBmin"   OBIT_float (1,1,1) min. beam gain for Squint correction [Jinc 0.05, Poly 0.01]
 * \li  "Factor"  OBIT_float (1,1,1) model multiplications factor (-1=>add) [def 1]
 * \li  "minFlux" OBIT_float (1,1,1) Minimum flux density model or pixel [def -1.0e20]
 * \li  "ModelType" OBIT_int (1,1,1) Model type (ObitSkyModelVMBeamType) [def OBIT_SkyModel_Comps]
 * \li  "Mode"     OBIT_int (1,1,1) Model mode (ObitSkyModelVMBeamMode) [def OBIT_SkyModel_Fastest]
 * \li  "MODPTFLX" OBIT_float (1,1,1) Point model flux in Jy, [def 0.0]
 * \li  "MODPTXOF" OBIT_float (1,1,1) Point model "x" offset in deg  [def 0.0]
 * \li  "MODPTYOF" OBIT_float (1,1,1) Point model "y" offset in deg  [def 0.0]
 * \li  "MODPTYPM" OBIT_float (4,1,1) Point other parameters  [def all 0.0]
 * \li  "CCVer" OBIT_int (?,1,1) CC table versions to use [def all 0 => highest]
 * \li  "BComp" OBIT_int (?,1,1) Start CC to use per table, 1-rel [def 1 ]
 * \li  "EComp" OBIT_int (?,1,1) Highest CC to use per table, 1-rel [def to end ]
 * \li  "UpdateInt" OBIT_float (1,1,1) Model update interval (min)  [def 1 min]
 * \li  "Threshold" OBIT_float (1,1,1) Threshold for high accuracy model [1.0e20]
 * \li  "BeamCorClean" OBIT_bool (1,1,1) Is this part of a CLEAN?
 * \li  "maxResid"  OBIT_float (1,1,1) Current maximum abs residual flux density
 *       if maxResid==0.0 then use accurate model for merged CC abs flux > Threshold and
 *       the gridded method for weaker summed components.
 *       [default or < 0 then determine from ImageMosaic]
 */

/*-------------- enumerations -------------------------------------*/

/*--------------Class definitions-------------------------------------*/
/** ObitSkyModelVMBeam Class structure. */
typedef struct {
#include "ObitSkyModelVMBeamDef.h"   /* this class definition */
} ObitSkyModelVMBeam;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitSkyModelVMBeam
 * returns a ObitSkyModelVMBeam*.
 * in = object to unreference
 */
#define ObitSkyModelVMBeamUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitSkyModelVMBeam.
 * returns a ObitSkyModelVMBeam*.
 * in = object to reference
 */
#define ObitSkyModelVMBeamRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitSkyModelVMBeamIsA(in) ObitIsA (in, ObitSkyModelVMBeamGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitSkyModelVMBeamClassInit (void);

/** Public: Default Constructor. */
ObitSkyModelVMBeam* newObitSkyModelVMBeam (gchar* name);

/** Public: Init SkyModel object from description in an ObitInfoList */
void ObitSkyModelVMBeamFromInfo (ObitSkyModel *out, gchar *prefix, 
				 ObitInfoList *inList, ObitErr *err);

/** Public: Create/initialize ObitSkyModelVMBeam structures */
ObitSkyModelVMBeam* ObitSkyModelVMBeamCreate (gchar* name, ObitImageMosaic* mosaic,
					      ObitUV *uvData, olong numAntType,
					      ObitImage **RXBeam,   ObitImage **LYBeam, 
					      ObitImage **RLBeam,   ObitImage **LRBeam, 
					      ObitImage **RXBeamPh, ObitImage **LYBeamPh, 
					      ObitImage **RLBeamPh, ObitImage **LRBeamPh, 
					      ofloat *Diams, ObitErr *err);

/** Public: Set Antenna Types */
void  ObitSkyModelVMBeamSetAnt (ObitSkyModelVMBeam* in, ObitUV *uvdata, ObitErr *err);

/** Public: Get Inputs. */
void  ObitSkyModelVMBeamGetInput (ObitSkyModel* inn, ObitErr *err);

/** Public: initialize ObitSkyModelVMBeam structures */
void ObitSkyModelVMBeamInitMod (ObitSkyModel* in,  ObitUV *uvdata, ObitErr *err);

/** Public: shutdown ObitSkyModel processes */
void ObitSkyModelVMBeamShutDownMod (ObitSkyModel* in,  ObitUV *uvdata, 
				    ObitErr *err);

/** Public: initialize model for pass through data */
void ObitSkyModelVMBeamInitModel (ObitSkyModel* in, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitSkyModelVMBeamGetClass (void);

/** Public: Copy (deep) constructor. */
ObitSkyModelVMBeam* 
ObitSkyModelVMBeamCopy  (ObitSkyModelVMBeam *in, 
			 ObitSkyModelVMBeam *out, ObitErr *err);

/** Public: Copy structure. */
void ObitSkyModelVMBeamClone (ObitSkyModelVMBeam *in, ObitSkyModelVMBeam *out, 
			      ObitErr *err);

/** Public: Routine to update model */
void ObitSkyModelVMBeamUpdateModel (ObitSkyModelVM *in, ofloat time, olong suba,
				    ObitUV *uvdata, olong ithread, ObitErr *err);

/** Private: Chose model type */
void  ObitSkyModelVMBeamChose (ObitSkyModel* in, ObitUV* uvdata);

/** Public: Extract information about underlying structures to ObitInfoList */
void ObitSkyModelVMBeamGetInfo (ObitSkyModel *in, gchar *prefix, 
				ObitInfoList *outList, 
				ObitErr *err);

/** Public: Jones correct and average SkyModel components */
void ObitSkyModelVMBeamJonesCorSum (olong numComp, olong Stokes, gboolean isCirc, 
				    ofloat *compFlux, ofloat *sinArr, ofloat *cosArr,
				    ObitMatx **Jones1, ObitMatx **Jones2, ofloat parAng2,
				    ObitMatx **workVis, ObitMatx *work1, ObitMatx *work2,
				    ObitMatx *sumVis);
/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitSkyModelVMBeamClassDef.h"
} ObitSkyModelVMBeamClassInfo; 

#endif /* OBITSKYMODELVMBEAM_H */ 
