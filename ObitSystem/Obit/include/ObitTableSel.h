/* $Id: ObitTableSel.h,v 1.3 2004/12/28 14:40:49 bcotton Exp $                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITTABLESEL_H 
#define OBITTABLESEL_H 
#include <glib.h>
#include "Obit.h"
#include "ObitTableDesc.h"
#include "ObitErr.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableSel.h
 * ObitTableSel Obit Table selector class definition.
 * This class is derived from the Obit class.
 * This contains the descriptions of data selection and calibration.
 *
 * \section ObitTableSelUsage Usage
 * Instances can be obtained using the #newObitTableSel constructor
 * the #ObitTableSelCopy copy constructor or a pointer duplicated using 
 * the #ObitTableSelRef function.
 * When an instance is no longer needed, use the #ObitTableSelUnref macro
 * to release it.
 */

/*------------------- Macroes ----------------------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitTableSel
 * returns a ObitTableSel* (NULL).
 * \li in = object to unreference.
 */
#define ObitTableSelUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitTableSel.
 * returns a ObitTableSel*.
 * in = object to reference
 */
#define ObitTableSelRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitTableSelIsA(in) ObitIsA (in, ObitTableSelGetClass())

/*--------------Class definitions-------------------------------------*/
/**
 * ObitTableSel Class structure.
 *
 * This class contains descriptions of interferometric visibility data.
 */  
typedef struct {
#include "ObitTableSelDef.h" /* actual definition */
} ObitTableSel;

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitTableSelClassInit (void);

/** Public: Constructor. */
ObitTableSel* newObitTableSel (gchar *name);

/** Public: Return class pointer. */
gconstpointer ObitTableSelGetClass (void);

/** Public: Copy TableSel */
ObitTableSel* ObitTableSelCopy (ObitTableSel* in, ObitTableSel* out,
			  ObitErr *err);

/** Public: How big a buffer is needed for a data transfer? */
olong ObitTableSelBufferSize (ObitTableDesc* desc, 
			      ObitTableSel* sel);

/** Public: Enforces defaults in inaxes, blc, trc */
void ObitTableSelDefault (ObitTableDesc* in, 
			  ObitTableSel* sel);

/** Public: Applies selection to a Descriptor */
void ObitTableSelSetDesc (ObitTableDesc* in, ObitTableSel* sel,
			  ObitTableDesc* out, ObitErr *err);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitTableSelClassDef.h" /* Actual definition */
} ObitTableSelClassInfo; 

#endif /* OBITTABLESEL_H */ 

