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
#ifndef OBITGBTVEGASSTATEINFO_H 
#define OBITGBTVEGASSTATEINFO_H 
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGBTVEGASStateInfo.h
 * ObitGBTVEGASStateInfo class definition.
 * This class is derived from the Obit class.
 *
 * This class contains information about the GBT observing VEGASState.
 * Frequency and polarization information are retrieved from GBT 
 * archive files.
 *
 * \section ObitGBTVEGASStateInfoUsage Usage
 * Instances can be obtained using the #newObitGBTVEGASStateInfo constructor,
 * the #ObitGBTVEGASStateInfoCopy constructor or a pointer duplicated using 
 * the #ObitGBTVEGASStateInfoRef macro.
 * When an instance is no longer needed, use the #ObitGBTVEGASStateInfoUnref 
 * macro to release it.
 */

/*---------------Class Structure---------------------------*/
/** ObitGBTVEGASStateInfo Class. */
typedef struct {
#include "ObitGBTVEGASStateInfoDef.h"   /* actual definition */
} ObitGBTVEGASStateInfo;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitGBTVEGASStateInfo
 * returns a ObitGBTVEGASStateInfo*.
 * in = object to unreference
 */
#define ObitGBTVEGASStateInfoUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitGBTVEGASStateInfo.
 * returns a ObitGBTVEGASStateInfo*.
 * in = object to reference
 */
#define ObitGBTVEGASStateInfoRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGBTVEGASStateInfoIsA(in) ObitIsA (in, ObitGBTVEGASStateInfoGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
 void ObitGBTVEGASStateInfoClassInit (void);

/** Public: Constructor. */
ObitGBTVEGASStateInfo* newObitGBTVEGASStateInfo (gchar* name);

/** Public: Constructor from values. */
ObitGBTVEGASStateInfo* 
newObitGBTVEGASStateInfoValue (gchar *name, gchar *DataRoot, gchar *scan, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitGBTVEGASStateInfoGetClass (void);

/** Public: Copy  constructor. */
ObitGBTVEGASStateInfo* 
ObitGBTVEGASStateInfoCopy  (ObitGBTVEGASStateInfo *in, ObitGBTVEGASStateInfo *out, ObitErr *err);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * and function pointers.
 */
typedef struct  {
#include "ObitGBTVEGASStateInfoClassDef.h" /* Actual definition */
} ObitGBTVEGASStateInfoClassInfo; 


#endif /* OBITGBTVEGASSTATEINFO_H */ 
