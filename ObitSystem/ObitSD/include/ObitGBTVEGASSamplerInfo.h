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
#ifndef OBITGBTVEGASSAMPLERINFO_H 
#define OBITGBTVEGASSAMPLERINFO_H 
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGBTVEGASSamplerInfo.h
 * ObitGBTVEGASSamplerInfo class definition.
 * This class is derived from the Obit class.
 *
 * This class contains information about the GBT observing VEGASSampler.
 * Frequency and polarization information are retrieved from GBT 
 * archive files.
 *
 * \section ObitGBTVEGASSamplerInfoUsage Usage
 * Instances can be obtained using the #newObitGBTVEGASSamplerInfo constructor,
 * the #ObitGBTVEGASSamplerInfoCopy constructor or a pointer duplicated using 
 * the #ObitGBTVEGASSamplerInfoRef macro.
 * When an instance is no longer needed, use the #ObitGBTVEGASSamplerInfoUnref 
 * macro to release it.
 */

/*---------------Class Structure---------------------------*/
/** ObitGBTVEGASSamplerInfo Class. */
typedef struct {
#include "ObitGBTVEGASSamplerInfoDef.h"   /* actual definition */
} ObitGBTVEGASSamplerInfo;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitGBTVEGASSamplerInfo
 * returns a ObitGBTVEGASSamplerInfo*.
 * in = object to unreference
 */
#define ObitGBTVEGASSamplerInfoUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitGBTVEGASSamplerInfo.
 * returns a ObitGBTVEGASSamplerInfo*.
 * in = object to reference
 */
#define ObitGBTVEGASSamplerInfoRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGBTVEGASSamplerInfoIsA(in) ObitIsA (in, ObitGBTVEGASSamplerInfoGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
 void ObitGBTVEGASSamplerInfoClassInit (void);

/** Public: Constructor. */
ObitGBTVEGASSamplerInfo* newObitGBTVEGASSamplerInfo (gchar* name);

/** Public: Constructor from values. */
ObitGBTVEGASSamplerInfo* 
newObitGBTVEGASSamplerInfoValue (gchar *name, gchar *DataRoot, gchar *scan, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitGBTVEGASSamplerInfoGetClass (void);

/** Public: Copy  constructor. */
ObitGBTVEGASSamplerInfo* 
ObitGBTVEGASSamplerInfoCopy  (ObitGBTVEGASSamplerInfo *in, ObitGBTVEGASSamplerInfo *out, ObitErr *err);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * and function pointers.
 */
typedef struct  {
#include "ObitGBTVEGASSamplerInfoClassDef.h" /* Actual definition */
} ObitGBTVEGASSamplerInfoClassInfo; 


#endif /* OBITGBTVEGASSAMPLERINFO_H */ 
