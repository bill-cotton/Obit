/* $Id$ */
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
#ifndef OBITGBTPORTINFO_H 
#define OBITGBTPORTINFO_H 
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGBTCCBPortInfo.h
 * ObitGBTCCBPortInfo class definition.
 * This class is derived from the Obit class.
 *
 * This class contains information about the GBT observing CCBPort.
 * Frequency and polarization information are retrieved from GBT 
 * archive files.
 *
 * \section ObitGBTCCBPortInfoUsage Usage
 * Instances can be obtained using the #newObitGBTCCBPortInfo constructor,
 * the #ObitGBTCCBPortInfoCopy constructor or a pointer duplicated using 
 * the #ObitGBTCCBPortInfoRef macro.
 * When an instance is no longer needed, use the #ObitGBTCCBPortInfoUnref 
 * macro to release it.
 */

/*---------------Class Structure---------------------------*/
/** ObitGBTCCBPortInfo Class. */
typedef struct {
#include "ObitGBTCCBPortInfoDef.h"   /* actual definition */
} ObitGBTCCBPortInfo;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitGBTCCBPortInfo
 * returns a ObitGBTCCBPortInfo*.
 * in = object to unreference
 */
#define ObitGBTCCBPortInfoUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitGBTCCBPortInfo.
 * returns a ObitGBTCCBPortInfo*.
 * in = object to reference
 */
#define ObitGBTCCBPortInfoRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGBTCCBPortInfoIsA(in) ObitIsA (in, ObitGBTCCBPortInfoGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
 void ObitGBTCCBPortInfoClassInit (void);

/** Public: Constructor. */
ObitGBTCCBPortInfo* newObitGBTCCBPortInfo (gchar* name);

/** Public: Constructor from values. */
ObitGBTCCBPortInfo* 
newObitGBTCCBPortInfoValue (gchar *name, olong disk, gchar *scan, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitGBTCCBPortInfoGetClass (void);

/** Public: Copy  constructor. */
ObitGBTCCBPortInfo* 
ObitGBTCCBPortInfoCopy  (ObitGBTCCBPortInfo *in, ObitGBTCCBPortInfo *out, ObitErr *err);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * and function pointers.
 */
typedef struct  {
#include "ObitGBTCCBPortInfoClassDef.h" /* Actual definition */
} ObitGBTCCBPortInfoClassInfo; 


#endif /* OBITGBTPORTINFO_H */ 
