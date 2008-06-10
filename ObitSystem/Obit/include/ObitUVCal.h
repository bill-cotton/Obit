/* $Id: ObitUVCal.h,v 1.5 2007/08/31 17:24:49 bcotton Exp $        */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITUVCAL_H 
#define OBITUVCAL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitUVDesc.h"
#include "ObitUVSel.h"
#include "ObitUVCalCalibrateDef.h"
#include "ObitUVCalBaselineDef.h"
#include "ObitUVCalFlagDef.h"
#include "ObitUVCalBandpassDef.h"
#include "ObitUVCalPolarizationDef.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVCal.h
 * ObitUVCal data calibration class definition.
 *
 * This class is derived from the #Obit class.
 *
 * This class is for the calibration, editing and selection of uv data.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitUVCal Class structure. */
typedef struct {
#include "ObitUVCalDef.h"   /* this class definition */
} ObitUVCal;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitUVCal
 * returns a ObitUVCal*.
 * in = object to unreference
 */
#define ObitUVCalUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitUVCal.
 * returns a ObitUVCal*.
 * in = object to reference
 */
#define ObitUVCalRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitUVCalIsA(in) ObitIsA (in, ObitUVCalGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitUVCalClassInit (void);

/** Public: Constructor. */
ObitUVCal* newObitUVCal (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitUVCalGetClass (void);

/** Public: Copy (deep) constructor. */
ObitUVCal* ObitUVCalCopy  (ObitUVCal *in, ObitUVCal *out, 
			   ObitErr *err);

/** Public: Copy (shallow) constructor. */
ObitUVCal* ObitUVCalClone (ObitUVCal *in, ObitUVCal *out);

/** Public: Initialize calibration */
void ObitUVCalStart (ObitUVCal *in, ObitUVSel *sel, ObitUVDesc *inDesc, 
		    ObitUVDesc *outDesc, ObitErr *err);
typedef void (*ObitUVCalStartFP) (ObitUVCal *in, ObitUVSel *sel, 
				 ObitUVDesc *inDesc, ObitUVDesc *outDesc, ObitErr *err);

/** Public: Apply Calibration */
gboolean  ObitUVCalApply (ObitUVCal *in, ofloat *visIn, ofloat *visOut, 
			  ObitErr *err);
typedef gboolean (*ObitUVCalApplyFP) (ObitUVCal *in, ofloat *visIn, 
				      ofloat *visOut, ObitErr *err);

/** Public: Shutdown */
void ObitUVCalShutdown (ObitUVCal *in, ObitErr *err);
typedef void (*ObitUVCalShutdownFP) (ObitUVCal *in, ObitErr *err);

/** Public: test if data wanted  */
gboolean ObitUVCalWant (ObitUVCal *in, ofloat time, olong ant1, olong ant2, 
			ofloat *RP, ofloat *visIn, ObitErr *err);
typedef gboolean (*ObitUVCalWantFP) (ObitUVCal *in, ofloat time, olong ant1, olong ant2, 
				     ofloat *RP, ofloat *visIn, ObitErr *err);

/** Public: Smooth data in frequency*/
void ObitUVCalSmooth (ObitUVCal *in, float time, olong ant1, olong ant2, 
		      ofloat *RP, ofloat *visIn, ObitErr *err);
typedef void (*ObitUVCalSmoothFP) (ObitUVCal *in, float time, olong ant1, olong ant2, 
				   ofloat *visIn, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitUVCalClassDef.h"
} ObitUVCalClassInfo; 

#endif /* OBITUVCAL_H */ 
