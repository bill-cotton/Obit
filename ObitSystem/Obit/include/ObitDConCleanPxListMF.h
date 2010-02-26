/* $Id: ObitDConCleanPxListMF.h 128 2009-09-23 14:48:29Z bill.cotton $ */
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
#ifndef OBITDCONCLEANPXLISTMF_H 
#define OBITDCONCLEANPXLISTMF_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitImageMosaic.h"
#include "ObitDConCleanWindow.h"
#include "ObitTableCC.h"
#include "ObitDConCleanPxList.h"
#include "ObitUV.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanPxListMF.h
 * ObitDConCleanPxListMF CLEAN image pixel list/CLEAN class.
 *
 * ObitDConCleanPxListMF  Visibility-based (Cotton-Schwab) Pixel List class 
 * for wideband imaging in the spectral imaging style
 * Fit S = S_0 exp(alpha log(nu/nu_0) + beta log(nu/nu_0)*2)
 * This class is derived from the #ObitPxList class.
 * 
 * \section ObitDConCleanPxListMFaccess Creators and Destructors
 * An ObitDConCleanPxListMF will usually be created using ObitDConCleanPxListMFCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitDConCleanPxListMF should always be made using the
 * #ObitDConCleanPxListMFRef function which updates the reference count in the object.
 * Then whenever freeing an ObitDConCleanPxListMF or changing a pointer, the function
 * #ObitDConCleanPxListMFUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*-------------- enumerations -------------------------------------*/

/*--------------Class definitions-------------------------------------*/
/** ObitDConCleanPxListMF Class structure. */
typedef struct {
#include "ObitDConCleanPxListMFDef.h"   /* this class definition */
} ObitDConCleanPxListMF;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitDConCleanPxListMF
 * returns a ObitDConCleanPxListMF*.
 * in = object to unreference
 */
#define ObitDConCleanPxListMFUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitDConCleanPxListMF.
 * returns a ObitDConCleanPxListMF*.
 * in = object to reference
 */
#define ObitDConCleanPxListMFRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitDConCleanPxListMFIsA(in) ObitIsA (in, ObitDConCleanPxListMFGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitDConCleanPxListMFClassInit (void);

/** Public: Default Constructor. */
ObitDConCleanPxListMF* newObitDConCleanPxListMF (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitDConCleanPxListMFGetClass (void);

/** Public: Copy (deep) constructor. */
ObitDConCleanPxListMF* 
ObitDConCleanPxListMFCopy  (ObitDConCleanPxListMF *in, ObitDConCleanPxListMF *out, 
			  ObitErr *err);

/** Public: Copy structure. */
void ObitDConCleanPxListMFClone (ObitDConCleanPxListMF *in, 
			       ObitDConCleanPxListMF *out, 
			       ObitErr *err);

/** Public: Create/initialize ObitDCon structures */
ObitDConCleanPxListMF* 
ObitDConCleanPxListMFCreate (gchar* name, ObitImageMosaic *mosaic, 
			     ObitUV *uvdata, olong maxPixel, ObitErr *err);

/** Public: Get Parameters. */
void ObitDConCleanPxListMFGetParms (ObitDConCleanPxList *in, ObitErr *err);

/** Public: Reset Clean. */
void ObitDConCleanPxListMFReset (ObitDConCleanPxList *in, ObitErr *err);

/** Public: Resize Arrrays. */
void ObitDConCleanPxListMFResize (ObitDConCleanPxList *in, olong maxPixel, 
				ObitErr *err);

/** Public: Update with new image and window. */
void ObitDConCleanPxListMFUpdate (ObitDConCleanPxList *in, 
				olong *fields, olong nSkip, 
				ofloat minFluxLoad,
				ofloat autoWinFlux,
				ObitDConCleanWindow *window, 
				ObitFArray **BeamPatch,
				ObitFArray **pixarray,
				ObitErr *err);

/** Public: Do minor cycle BGC CLEANing. */
gboolean ObitDConCleanPxListMFCLEAN (ObitDConCleanPxList *in, ObitErr *err);

/** Public: Do SDI CLEANing. */
gboolean ObitDConCleanPxListMFSDI (ObitDConCleanPxList *in, ObitErr *err);

/** Public: Get results of CLEAN */
olong ObitDConCleanPxListMFResult (ObitDConCleanPxList *in, olong *ncomp,
				 ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitDConCleanPxListMFClassDef.h"
} ObitDConCleanPxListMFClassInfo; 

#endif /*  OBITDCONCLEANBMLISTMF_H  */ 
