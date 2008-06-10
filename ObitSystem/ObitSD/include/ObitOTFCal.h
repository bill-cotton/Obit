/* $Id$      */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITOTFCAL_H 
#define OBITOTFCAL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitOTFDesc.h"
#include "ObitOTFSel.h"
#include "ObitOTFCalFlagDef.h"
#include "ObitOTFArrayGeom.h"
#include "ObitTableOTFCal.h"
#include "ObitTableOTFSoln.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFCal.h
 * ObitOTFCal "On the Fly" calibration class definition.
 * This class is derived from the Obit class.
 *
 * This class is for the calibration, editing and selection of OTF data.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitOTFCal Class structure. */
typedef struct {
#include "ObitOTFCalDef.h"   /* this class definition */
} ObitOTFCal;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitOTFCal
 * returns a ObitOTFCal*.
 * in = object to unreference
 */
#define ObitOTFCalUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitOTFCal.
 * returns a ObitOTFCal*.
 * in = object to reference
 */
#define ObitOTFCalRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitOTFCalIsA(in) ObitIsA (in, ObitOTFCalGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitOTFCalClassInit (void);

/** Public: Constructor. */
ObitOTFCal* newObitOTFCal (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitOTFCalGetClass (void);

/** Public: Copy (deep) constructor. */
ObitOTFCal* ObitOTFCalCopy  (ObitOTFCal *in, ObitOTFCal *out, 
			   ObitErr *err);

/** Public: Copy (shallow) constructor. */
ObitOTFCal* ObitOTFCalClone (ObitOTFCal *in, ObitOTFCal *out);

/** Public: Initialize calibration */
void ObitOTFCalStart (ObitOTFCal *in, ObitOTFSel *sel, ObitOTFDesc *inDesc, 
		      ObitOTFArrayGeom *geom, ObitOTFDesc *outDesc, ObitErr *err);
typedef void (*ObitOTFCalStartFP) (ObitOTFCal *in, ObitOTFSel *sel, 
				   ObitOTFArrayGeom *geom, ObitOTFDesc *inDesc, ObitOTFDesc *outDesc, 
				   ObitErr *err);

/** Public: Apply Calibration */
gboolean  ObitOTFCalApply (ObitOTFCal *in, ofloat *recIn, ofloat *recOut, ObitErr *err);
typedef gboolean (*ObitOTFCalApplyFP) (ObitOTFCal *in, ofloat *recIn, ofloat *recOut, ObitErr *err);

/** Public: Shutdown */
ObitOTFCal* ObitOTFCalShutdown (ObitOTFCal *in, ObitErr *err);
typedef void (*ObitOTFCalShutdownFP) (ObitOTFCal *in, ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitOTFCalClassDef.h"
} ObitOTFCalClassInfo; 

#endif /* OBITOTFCAL_H */ 
