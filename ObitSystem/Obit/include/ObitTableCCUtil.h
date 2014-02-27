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
#ifndef OBITTABLECCUTIL_H 
#define OBITTABLECCUTIL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitTableCC.h"
#include "ObitImage.h"
#include "ObitImageMF.h"
#include "ObitImageWB.h"
#include "ObitImageDesc.h"
#include "ObitFArray.h"
#include "ObitData.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableCCUtil.h
 * ObitTableCC class utility routine definition.
 */

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitCCCompType
 * enum for CC table Component model type
 */
enum obitCCCompType {
  /** Point */
  OBIT_CC_PointMod,
  /** Gaussian on sky */
  OBIT_CC_GaussMod,  
  /** Convolved Gaussian */
  OBIT_CC_CGaussMod,  
  /** Uniform sphere */
  OBIT_CC_USphereMod, 
  /** Unknown */
  OBIT_CC_Unknown, 
  /** Point + spectrum */
  OBIT_CC_PointModSpec = OBIT_CC_PointMod+10,
  /** Gaussian on sky + spectrum */
  OBIT_CC_GaussModSpec = OBIT_CC_GaussMod+10,  
  /** Convolved Gaussian + spectrum */
  OBIT_CC_CGaussModSpec = OBIT_CC_CGaussMod+10,  
  /** Uniform sphere + spectrum */
  OBIT_CC_USphereModSpec = OBIT_CC_USphereMod+10,
  /** Point + tabulated spectrum */
  OBIT_CC_PointModTSpec = OBIT_CC_PointMod+20,
  /** Gaussian on sky + tabulated spectrum */
  OBIT_CC_GaussModTSpec = OBIT_CC_GaussMod+20,  
  /** Convolved Gaussian + tabulated spectrum */
  OBIT_CC_CGaussModTSpec = OBIT_CC_CGaussMod+20,  
  /** Uniform sphere + tabulated spectrum */
  OBIT_CC_USphereModTSpec = OBIT_CC_USphereMod+20 
}; /* end enum obitCCCompType */
/** typedef for enum for ObitCCCompType. */
typedef enum obitCCCompType ObitCCCompType;

/*---------------Public functions---------------------------*/
/** Public: grid components onto a grid */
ObitIOCode ObitTableCCUtilGrid (ObitTableCC *in, olong OverSample, 
				olong *first, olong *last, gboolean noNeg,
				ofloat factor, ofloat minFlux, ofloat maxFlux,
				ObitImageDesc *desc, ObitFArray **grid, 
				ofloat gparm[3], olong *ncomps,
				ObitErr *err);

/** Public: grid spectral components onto a grid */
ObitIOCode ObitTableCCUtilGridSpect (ObitTableCC *in, olong OverSample, olong iterm,
				olong *first, olong *last, gboolean noNeg,
				ofloat factor, ofloat minFlux, ofloat maxFlux,
				ObitImageDesc *desc, ObitFArray **grid, 
				ofloat gparm[3], olong *ncomps,
				ObitErr *err);

/** Public: return list of CC from one image overlapping another */
ObitFArray* 
ObitTableCCUtilCrossList (ObitTableCC *inCC, ObitImageDesc *inDesc,  
			  ObitImageDesc *outDesc, ofloat gparm[3], 
			  olong *ncomps, ObitErr *err);

/** Public: return list of spectral components from one image overlapping another */
ObitFArray* 
ObitTableCCUtilCrossListSpec (ObitTableCC *inCC, ObitImageDesc *inDesc,  
			      ObitImageDesc *outDesc, ofloat gparm[3], 
			      olong *ncomps, olong iterm, ObitErr *err);

/** Public: return Table of CC from one image overlapping another */
ObitTableCC* 
ObitTableCCUtilCrossTable (ObitTableCC *inCC, ObitImageDesc *inDesc,  
			   ObitImage *outIm, olong *ncomps, 
			   ObitErr *err);

/** Merge elements of an ObitTableCC */
ObitIOCode ObitTableCCUtilMerge (ObitTableCC *in, ObitTableCC *out, 
				 ObitErr *err);

/** Merge elements of an ObitTableCC with selection */
ObitFArray* ObitTableCCUtilMergeSel (ObitTableCC *in, olong startComp, 
				     olong endComp, ofloat *parms,
				     ObitErr *err);

/** Merge spectral elements of an ObitTableCC with selection */
ObitFArray* ObitTableCCUtilMergeSelSpec (ObitTableCC *in, olong startComp, 
					 olong endComp, ofloat *parms,
					 ObitErr *err);

/** Merge selected elements of an ObitTableCC to a new table */
ObitTableCC* 
ObitTableCCUtilMergeSel2Tab (ObitImage *image, olong inCCver, olong *outCCver,
			     olong startComp, olong endComp, 
			     ofloat range[2], ObitErr *err);

/** Scale the flux densities of entries in a CC table */
void ObitTableCCUtilScale (ObitTableCC *in, olong startComp, 
			   olong endComp, ofloat scale, ObitErr *err);

/** Append CLEAN components from one table to another */
void ObitTableCCUtilAppend  (ObitTableCC *inCC, ObitTableCC *outCC, 
			     olong startComp, olong endComp, ObitErr *err);

/** Append CLEAN components from one table to another with position shift */
void ObitTableCCUtilAppendShift (ObitTableCC *inCC, ObitTableCC *outCC, 
				 ObitUVDesc *uvDesc, ObitImageDesc *imDesc, 
				 olong startComp, olong endComp, ObitErr *err);

/** Filter weak, isolated components */
gboolean ObitTableCCUtilFiltCC (ObitTableCC *CCTab, ofloat radius, ofloat minFlux, 
				ObitErr* err);

/** Get Clean component type */
ObitCCCompType ObitTableCCUtilGetType (ObitData *data, olong ver, ObitErr* err);

/** Convert TSpec to Spec model type */
void ObitTableCCUtilT2Spec  (ObitImage *inImage, ObitImageWB *outImage, 
			     olong nTerm, olong *inCCVer, olong *outCCVer,
			     olong startComp, olong endComp, ObitErr *err);

/* routine to force average TSpectra to a given spectrum  */
void ObitTableCCUtilFixTSpec (ObitImage *inImage, olong *inCCVer, 
			      odouble refFreq, olong nterm, ofloat *terms,
			      olong startCC, olong endCC, ObitErr *err);

/* routine to combine simple CC Tables into a TSpec table  */
void ObitTableCCUtilCombTSpec (ObitImage *inImage, olong inCCVer, olong nCCVer,
			       olong outCCVer, olong *bcopy, olong *bcomp, olong *ecomp, 
			       gboolean doGaus, ObitErr *err);
#endif /* OBITTABLECCUTIL_H */ 
