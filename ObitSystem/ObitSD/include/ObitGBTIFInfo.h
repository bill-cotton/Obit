/* $Id$                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2013                                          */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITGBTIFINFO_H 
#define OBITGBTIFINFO_H 
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGBTIFInfo.h
 * ObitGBTIFInfo class definition.
 * This class is derived from the Obit class.
 *
 * This class contains information about a GBT IF setup.
 * Frequency and polarization information are retrieved from GBT 
 * archive files.
 *
 * \section ObitGBTIFInfoUsage Usage
 * Instances can be obtained using the #newObitGBTIFInfo constructor,
 * the #ObitGBTIFInfoCopy constructor or a pointer duplicated using 
 * the #ObitGBTIFInfoRef macro.
 * When an instance is no longer needed, use the #ObitGBTIFInfoUnref 
 * macro to release it.
 */

/*---------------Class Structure---------------------------*/
/** ObitGBTIFInfo Class. */
typedef struct {
#include "ObitGBTIFInfoDef.h"   /* actual definition */
} ObitGBTIFInfo;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitGBTIFInfo
 * returns a ObitGBTIFInfo*.
 * in = object to unreference
 */
#define ObitGBTIFInfoUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitGBTIFInfo.
 * returns a ObitGBTIFInfo*.
 * in = object to reference
 */
#define ObitGBTIFInfoRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGBTIFInfoIsA(in) ObitIsA (in, ObitGBTIFInfoGetClass())

/*---------------Public functions---------------------------*/
/**  Public: Class initializer. */
void ObitGBTIFInfoClassInit (void);

/** Public: Constructor. */
ObitGBTIFInfo* newObitGBTIFInfo (gchar* name);

/** Public: Constructor from values. */
ObitGBTIFInfo* 
newObitGBTIFInfoValue (gchar *name, gchar *backend, olong disk, gchar *scan, ObitErr *err);
ObitGBTIFInfo* 
newObitGBTIFInfoValueRoot (gchar *name, gchar *backend, gchar *DataRoot, gchar *scan, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitGBTIFInfoGetClass (void);

/** Public: Copy  constructor. */
ObitGBTIFInfo* 
ObitGBTIFInfoCopy  (ObitGBTIFInfo *in, ObitGBTIFInfo *out, ObitErr *err);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * and function pointers.
 */
typedef struct  {
#include "ObitGBTIFInfoClassDef.h" /* Actual definition */
} ObitGBTIFInfoClassInfo; 


#endif /* OBITGBTIFINFO_H */ 
