/* $Id$                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004                                               */
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
 * \file ObitGBTDCRStateInfo.h
 * ObitGBTDCRStateInfo class definition.
 * This class is derived from the Obit class.
 *
 * This class contains information about the GBT observing DCRState.
 * Frequency and polarization information are retrieved from GBT 
 * archive files.
 *
 * \section ObitGBTDCRStateInfoUsage Usage
 * Instances can be obtained using the #newObitGBTDCRStateInfo constructor,
 * the #ObitGBTDCRStateInfoCopy constructor or a pointer duplicated using 
 * the #ObitGBTDCRStateInfoRef macro.
 * When an instance is no longer needed, use the #ObitGBTDCRStateInfoUnref 
 * macro to release it.
 */

/*---------------Class Structure---------------------------*/
/** ObitGBTDCRStateInfo Class. */
typedef struct {
#include "ObitGBTDCRStateInfoDef.h"   /* actual definition */
} ObitGBTDCRStateInfo;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitGBTDCRStateInfo
 * returns a ObitGBTDCRStateInfo*.
 * in = object to unreference
 */
#define ObitGBTDCRStateInfoUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitGBTDCRStateInfo.
 * returns a ObitGBTDCRStateInfo*.
 * in = object to reference
 */
#define ObitGBTDCRStateInfoRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGBTDCRStateInfoIsA(in) ObitIsA (in, ObitGBTDCRStateInfoGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
 void ObitGBTDCRStateInfoClassInit (void);

/** Public: Constructor. */
ObitGBTDCRStateInfo* newObitGBTDCRStateInfo (gchar* name);

/** Public: Constructor from values. */
ObitGBTDCRStateInfo* 
newObitGBTDCRStateInfoValue (gchar *name, olong disk, gchar *scan, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitGBTDCRStateInfoGetClass (void);

/** Public: Copy  constructor. */
ObitGBTDCRStateInfo* 
ObitGBTDCRStateInfoCopy  (ObitGBTDCRStateInfo *in, ObitGBTDCRStateInfo *out, ObitErr *err);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * and function pointers.
 */
typedef struct  {
#include "ObitGBTDCRStateInfoClassDef.h" /* Actual definition */
} ObitGBTDCRStateInfoClassInfo; 


#endif /* OBITGBTSTATEINFO_H */ 
