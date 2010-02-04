/* $Id$     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2010                                          */
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
#include "ObitImageDesc.h"
#include "ObitFArray.h"

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
  OBIT_CC_USphereModSpec = OBIT_CC_USphereMod+10 
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

/** Filter weak, isolated components */
gboolean ObitTableCCUtilFiltCC (ObitTableCC *CCTab, ofloat radius, ofloat minFlux, 
				ObitErr* err);
#endif /* OBITTABLECCUTIL_H */ 
