/* $Id: ObitGBTBeamOffInfo.h,v 1.3 2004/12/28 13:52:44 bcotton Exp $*/
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
#ifndef OBITGBTBEAMOFFINFO_H 
#define OBITGBTBEAMOFFINFO_H 
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGBTBeamOffInfo.h
 * ObitGBTBeamOffInfo class definition.
 * This class is derived from the Obit class.
 *
 * This class contains information about GBT Beam Offsets
 * Frequency and polarization information are retrieved from GBT 
 * archive files.
 *
 * \section ObitGBTBeamOffInfoUsage Usage
 * Instances can be obtained using the #newObitGBTBeamOffInfo constructor,
 * the #ObitGBTBeamOffInfoCopy constructor or a pointer duplicated using 
 * the #ObitGBTBeamOffInfoRef macro.
 * When an instance is no longer needed, use the #ObitGBTBeamOffInfoUnref 
 * macro to release it.
 */

/*---------------Class Structure---------------------------*/
/** ObitGBTBeamOffInfo Class. */
typedef struct {
#include "ObitGBTBeamOffInfoDef.h"   /* actual definition */
} ObitGBTBeamOffInfo;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitGBTBeamOffInfo
 * returns a ObitGBTBeamOffInfo*.
 * in = object to unreference
 */
#define ObitGBTBeamOffInfoUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitGBTBeamOffInfo.
 * returns a ObitGBTBeamOffInfo*.
 * in = object to reference
 */
#define ObitGBTBeamOffInfoRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGBTBeamOffInfoIsA(in) ObitIsA (in, ObitGBTBeamOffInfoGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitGBTBeamOffInfoClassInit (void);

/** Public: Constructor. */
ObitGBTBeamOffInfo* newObitGBTBeamOffInfo (gchar* name);

/** Public: Constructor from values. */
ObitGBTBeamOffInfo* 
newObitGBTBeamOffInfoValue (gchar *name, olong disk, gchar *scan, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitGBTBeamOffInfoGetClass (void);

/** Public: Copy  constructor. */
ObitGBTBeamOffInfo* 
ObitGBTBeamOffInfoCopy  (ObitGBTBeamOffInfo *in, ObitGBTBeamOffInfo *out, ObitErr *err);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * and function pointers.
 */
typedef struct  {
#include "ObitGBTBeamOffInfoClassDef.h" /* Actual definition */
} ObitGBTBeamOffInfoClassInfo; 


#endif /* OBITGBTBEAMOFFINFO_H */ 
