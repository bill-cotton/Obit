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
#ifndef OBITGBTVEGASACT_STATEINFO_H 
#define OBITGBTVEGASACT_STATEINFO_H 
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGBTVEGASAct\_StateInfo.h
 * ObitGBTVEGASAct_StateInfo class definition.
 * This class is derived from the Obit class.
 *
 * This class contains information about the GBT observing VEGASAct\_State.
 * Frequency and polarization information are retrieved from GBT 
 * archive files.
 *
 * \section ObitGBTVEGASAct\_StateInfoUsage Usage
 * Instances can be obtained using the #newObitGBTVEGASAct_StateInfo constructor,
 * the #ObitGBTVEGASAct_StateInfoCopy constructor or a pointer duplicated using 
 * the #ObitGBTVEGASAct_StateInfoRef macro.
 * When an instance is no longer needed, use the #ObitGBTVEGASAct_StateInfoUnref 
 * macro to release it.
 */

/*---------------Class Structure---------------------------*/
/** ObitGBTVEGASAct_StateInfo Class. */
typedef struct {
#include "ObitGBTVEGASAct_StateInfoDef.h"   /* actual definition */
} ObitGBTVEGASAct_StateInfo;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitGBTVEGASAct\_StateInfo
 * returns a ObitGBTVEGASAct_StateInfo*.
 * in = object to unreference
 */
#define ObitGBTVEGASAct_StateInfoUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitGBTVEGASAct\_StateInfo.
 * returns a ObitGBTVEGASAct_StateInfo*.
 * in = object to reference
 */
#define ObitGBTVEGASAct_StateInfoRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGBTVEGASAct_StateInfoIsA(in) ObitIsA (in, ObitGBTVEGASAct_StateInfoGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
 void ObitGBTVEGASAct_StateInfoClassInit (void);

/** Public: Constructor. */
ObitGBTVEGASAct_StateInfo* newObitGBTVEGASAct_StateInfo (gchar* name);

/** Public: Constructor from values. */
ObitGBTVEGASAct_StateInfo* 
newObitGBTVEGASAct_StateInfoValue (gchar *name, gchar *DataRoot, 
				   gchar *scan, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitGBTVEGASAct_StateInfoGetClass (void);

/** Public: Copy  constructor. */
ObitGBTVEGASAct_StateInfo* 
ObitGBTVEGASAct_StateInfoCopy  (ObitGBTVEGASAct_StateInfo *in, 
				ObitGBTVEGASAct_StateInfo *out, ObitErr *err);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * and function pointers.
 */
typedef struct  {
#include "ObitGBTVEGASAct_StateInfoClassDef.h" /* Actual definition */
} ObitGBTVEGASAct_StateInfoClassInfo; 


#endif /* OBITGBTVEGASACT_STATEINFO_H */ 
