/* $Id: ObitGBTCCBStateInfo.h,v 1.1 2006/03/22 18:47:11 bcotton Exp $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006                                               */
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
#ifndef OBITGBTSTATEINFO_H 
#define OBITGBTSTATEINFO_H 
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGBTCCBStateInfo.h
 * ObitGBTCCBStateInfo class definition.
 * This class is derived from the Obit class.
 *
 * This class contains information about the GBT observing CCBState.
 * Frequency and polarization information are retrieved from GBT 
 * archive files.
 *
 * \section ObitGBTCCBStateInfoUsage Usage
 * Instances can be obtained using the #newObitGBTCCBStateInfo constructor,
 * the #ObitGBTCCBStateInfoCopy constructor or a pointer duplicated using 
 * the #ObitGBTCCBStateInfoRef macro.
 * When an instance is no longer needed, use the #ObitGBTCCBStateInfoUnref 
 * macro to release it.
 */

/*---------------Class Structure---------------------------*/
/** ObitGBTCCBStateInfo Class. */
typedef struct {
#include "ObitGBTCCBStateInfoDef.h"   /* actual definition */
} ObitGBTCCBStateInfo;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitGBTCCBStateInfo
 * returns a ObitGBTCCBStateInfo*.
 * in = object to unreference
 */
#define ObitGBTCCBStateInfoUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitGBTCCBStateInfo.
 * returns a ObitGBTCCBStateInfo*.
 * in = object to reference
 */
#define ObitGBTCCBStateInfoRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGBTCCBStateInfoIsA(in) ObitIsA (in, ObitGBTCCBStateInfoGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
 void ObitGBTCCBStateInfoClassInit (void);

/** Public: Constructor. */
ObitGBTCCBStateInfo* newObitGBTCCBStateInfo (gchar* name);

/** Public: Constructor from values. */
ObitGBTCCBStateInfo* 
newObitGBTCCBStateInfoValue (gchar *name, olong disk, gchar *scan, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitGBTCCBStateInfoGetClass (void);

/** Public: Copy  constructor. */
ObitGBTCCBStateInfo* 
ObitGBTCCBStateInfoCopy  (ObitGBTCCBStateInfo *in, ObitGBTCCBStateInfo *out, ObitErr *err);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * and function pointers.
 */
typedef struct  {
#include "ObitGBTCCBStateInfoClassDef.h" /* Actual definition */
} ObitGBTCCBStateInfoClassInfo; 


#endif /* OBITGBTSTATEINFO_H */ 
