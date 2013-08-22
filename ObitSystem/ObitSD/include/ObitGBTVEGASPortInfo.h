/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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
#ifndef OBITGBTVEGASPORTINFO_H 
#define OBITGBTVEGASPORTINFO_H 
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGBTVEGASPortInfo.h
 * ObitGBTVEGASPortInfo class definition.
 * This class is derived from the Obit class.
 *
 * This class contains information about the GBT observing VEGASPort.
 * Frequency and polarization information are retrieved from GBT 
 * archive files.
 *
 * \section ObitGBTVEGASPortInfoUsage Usage
 * Instances can be obtained using the #newObitGBTVEGASPortInfo constructor,
 * the #ObitGBTVEGASPortInfoCopy constructor or a pointer duplicated using 
 * the #ObitGBTVEGASPortInfoRef macro.
 * When an instance is no longer needed, use the #ObitGBTVEGASPortInfoUnref 
 * macro to release it.
 */

/*---------------Class Structure---------------------------*/
/** ObitGBTVEGASPortInfo Class. */
typedef struct {
#include "ObitGBTVEGASPortInfoDef.h"   /* actual definition */
} ObitGBTVEGASPortInfo;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitGBTVEGASPortInfo
 * returns a ObitGBTVEGASPortInfo*.
 * in = object to unreference
 */
#define ObitGBTVEGASPortInfoUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitGBTVEGASPortInfo.
 * returns a ObitGBTVEGASPortInfo*.
 * in = object to reference
 */
#define ObitGBTVEGASPortInfoRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGBTVEGASPortInfoIsA(in) ObitIsA (in, ObitGBTVEGASPortInfoGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
 void ObitGBTVEGASPortInfoClassInit (void);

/** Public: Constructor. */
ObitGBTVEGASPortInfo* newObitGBTVEGASPortInfo (gchar* name);

/** Public: Constructor from values. */
ObitGBTVEGASPortInfo* 
newObitGBTVEGASPortInfoValue (gchar *name, gchar *DataRoot, gchar *scan, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitGBTVEGASPortInfoGetClass (void);

/** Public: Copy  constructor. */
ObitGBTVEGASPortInfo* 
ObitGBTVEGASPortInfoCopy  (ObitGBTVEGASPortInfo *in, ObitGBTVEGASPortInfo *out, ObitErr *err);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * and function pointers.
 */
typedef struct  {
#include "ObitGBTVEGASPortInfoClassDef.h" /* Actual definition */
} ObitGBTVEGASPortInfoClassInfo; 


#endif /* OBITGBTVEGASPORTINFO_H */ 
