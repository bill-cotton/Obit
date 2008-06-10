/* $Id: ObitSource.h,v 1.5 2007/08/31 17:24:48 bcotton Exp $    */
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
#ifndef OBITSOURCE_H 
#define OBITSOURCE_H 
#include "Obit.h"
#include "ObitErr.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSource.h
 * ObitSource class definition.
 *
 * This class is derived from the #Obit class.
 *
 * This class contains information about a given source.
 *
 * \section ObitSourceUsage Usage
 * Instances can be obtained using the #newObitSource constructor,
 * the #ObitSourceCopy constructor or a pointer duplicated using 
 * the #ObitSourceRef macro.
 * When an instance is no longer needed, use the #ObitSourceUnref 
 * macro to release it.
 */

/*---------------Class Structure---------------------------*/
/** ObitSource Class. */
typedef struct {
#include "ObitSourceDef.h"   /* actual definition */
} ObitSource;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitSource
 * returns a ObitSource*.
 * in = object to unreference
 */
#define ObitSourceUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitSource.
 * returns a ObitSource*.
 * in = object to reference
 */
#define ObitSourceRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitSourceIsA(in) ObitIsA (in, ObitSourceGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitSourceClassInit (void);

/** Public: Constructor. */
ObitSource* newObitSource (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitSourceGetClass (void);

/** Public: Copy  constructor. */
ObitSource* 
ObitSourceCopy  (ObitSource *in, ObitSource *out, ObitErr *err);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * and function pointers.
 */
typedef struct  {
#include "ObitSourceClassDef.h" /* Actual definition */
} ObitSourceClassInfo; 


#endif /* OBITSOURCE_H */ 
