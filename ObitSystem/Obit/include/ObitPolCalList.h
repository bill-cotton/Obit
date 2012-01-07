/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012                                               */
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
#ifndef OBITPOLCALLIST_H 
#define OBITPOLCALLIST_H 
#include "Obit.h"
#include "ObitErr.h"
#include "ObitAntenna.h"
#include "ObitAntennaList.h"
#include "ObitSource.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPolCalList.h
 * ObitPolCalList class definition.
 *
 * This class is derived from the #Obit class.
 *
 * This class manages lists of polarization calibration parameters.
 *
 * \section ObitPolCalListUsage Usage
 * Instances can be obtained using the #newObitPolCalList constructor,
 * the #ObitPolCalListCopy constructor or a pointer duplicated using 
 * the #ObitPolCalListRef macro.
 * When an instance is no longer needed, use the #ObitPolCalListUnref 
 * macro to release it.
 */

/*---------------Class Structure---------------------------*/
/** ObitPolCalList Class. */
typedef struct {
#include "ObitPolCalListDef.h"   /* actual definition */
} ObitPolCalList;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitPolCalList
 * returns a ObitPolCalList*.
 * in = object to unreference
 */
#define ObitPolCalListUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitPolCalList.
 * returns a ObitPolCalList*.
 * in = object to reference
 */
#define ObitPolCalListRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitPolCalListIsA(in) ObitIsA (in, ObitPolCalListGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
 void ObitPolCalListClassInit (void);

/** Public: Constructor. */
ObitPolCalList* newObitPolCalList (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitPolCalListGetClass (void);

/** Public: Determine polarization calibration type. */
ObitUVPolCalType ObitPolCalListGetPolType (gchar* type);

/** Public: Create from value */
ObitPolCalList* ObitPolCalListCreate (gchar* name, Obit *PDTab, 
				      ObitErr *err);


/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * and function pointers.
 */
typedef struct  {
#include "ObitPolCalListClassDef.h" /* Actual definition */
} ObitPolCalListClassInfo; 


#endif /* OBITPOLCALLIST_H */ 
