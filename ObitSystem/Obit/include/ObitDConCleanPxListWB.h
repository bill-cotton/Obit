/* $Id: ObitDConCleanPxListWB.h 128 2009-09-23 14:48:29Z bill.cotton $ */
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
#ifndef OBITDCONCLEANPXLISTWB_H 
#define OBITDCONCLEANPXLISTWB_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitImageMosaic.h"
#include "ObitDConCleanWindow.h"
#include "ObitTableCC.h"
#include "ObitDConCleanPxList.h"
#include "ObitUV.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanPxListWB.h
 * ObitDConCleanPxListWB CLEAN image pixel list/CLEAN class.
 *
 * ObitDConCleanPxListWB  Visibility-based (Cotton-Schwab) Pixel List class 
 * for wideband imaging in the style of Sault & Wieringa 1994 A&AS 108, 585
 * Fit S = S_0 exp(alpha log(nu/nu_0) + beta log(nu/nu_0)*2)
 * using approximation:
 * S = S0(1+ alpha*log(nu/nu_0) + beta * log(nu/nu_0)^2)
 * This class is derived from the #ObitPxList class.
 * 
 * \section ObitDConCleanPxListWBaccess Creators and Destructors
 * An ObitDConCleanPxListWB will usually be created using ObitDConCleanPxListWBCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitDConCleanPxListWB should always be made using the
 * #ObitDConCleanPxListWBRef function which updates the reference count in the object.
 * Then whenever freeing an ObitDConCleanPxListWB or changing a pointer, the function
 * #ObitDConCleanPxListWBUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*-------------- enumerations -------------------------------------*/

/*--------------Class definitions-------------------------------------*/
/** ObitDConCleanPxListWB Class structure. */
typedef struct {
#include "ObitDConCleanPxListWBDef.h"   /* this class definition */
} ObitDConCleanPxListWB;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitDConCleanPxListWB
 * returns a ObitDConCleanPxListWB*.
 * in = object to unreference
 */
#define ObitDConCleanPxListWBUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitDConCleanPxListWB.
 * returns a ObitDConCleanPxListWB*.
 * in = object to reference
 */
#define ObitDConCleanPxListWBRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitDConCleanPxListWBIsA(in) ObitIsA (in, ObitDConCleanPxListWBGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitDConCleanPxListWBClassInit (void);

/** Public: Default Constructor. */
ObitDConCleanPxListWB* newObitDConCleanPxListWB (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitDConCleanPxListWBGetClass (void);

/** Public: Copy (deep) constructor. */
ObitDConCleanPxListWB* 
ObitDConCleanPxListWBCopy  (ObitDConCleanPxListWB *in, ObitDConCleanPxListWB *out, 
			  ObitErr *err);

/** Public: Copy structure. */
void ObitDConCleanPxListWBClone (ObitDConCleanPxListWB *in, 
			       ObitDConCleanPxListWB *out, 
			       ObitErr *err);

/** Public: Create/initialize ObitDCon structures */
ObitDConCleanPxListWB* 
ObitDConCleanPxListWBCreate (gchar* name, ObitImageMosaic *mosaic, 
			     ObitUV *uvdata, olong maxPixel, ObitErr *err);

/** Public: Get Parameters. */
void ObitDConCleanPxListWBGetParms (ObitDConCleanPxList *in, ObitErr *err);

/** Public: Reset Clean. */
void ObitDConCleanPxListWBReset (ObitDConCleanPxList *in, ObitErr *err);

/** Public: Resize Arrrays. */
void ObitDConCleanPxListWBResize (ObitDConCleanPxList *in, olong maxPixel, 
				ObitErr *err);

/** Public: Update with new image and window. */
void ObitDConCleanPxListWBUpdate (ObitDConCleanPxList *in, 
				olong *fields, olong nSkip, 
				ofloat minFluxLoad,
				ofloat autoWinFlux,
				ObitDConCleanWindow *window, 
				ObitFArray **BeamPatch,
				ObitFArray **pixarray,
				ObitErr *err);

/** Public: Do minor cycle BGC CLEANing. */
gboolean ObitDConCleanPxListWBCLEAN (ObitDConCleanPxList *in, ObitErr *err);

/** Public: Do SDI CLEANing. */
gboolean ObitDConCleanPxListWBSDI (ObitDConCleanPxList *in, ObitErr *err);

/** Public: Get results of CLEAN */
olong ObitDConCleanPxListWBResult (ObitDConCleanPxList *in, olong *ncomp,
				 ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitDConCleanPxListWBClassDef.h"
} ObitDConCleanPxListWBClassInfo; 

#endif /*  OBITDCONCLEANBMLISTWB_H  */ 
