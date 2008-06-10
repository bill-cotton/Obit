/* $Id: ObitGBTDCROTF.h,v 1.4 2007/09/11 12:50:20 bcotton Exp $     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                          */
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
#ifndef OBITGBTOTF_H 
#define OBITGBTOTF_H 
#include "ObitOTF.h"
#include "Obit.h"
#include "ObitErr.h"
#include "ObitIOOTFFITS.h"
#include "ObitFITS.h"
#include "ObitTableOTFTarget.h"
#include "ObitTableOTFIndex.h"
#include "ObitTableGBTANTPOSGR.h"
#include "ObitTableGBTANTPOSPF.h"
#include "ObitTableGBTDCRDATA.h"
#include "ObitGBTIFInfo.h"
#include "ObitGBTDCRStateInfo.h"
#include "ObitGBTBeamOffInfo.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGBTDCROTF.h
 *
 * ObitGBTDCROTF class definition.
 * This class is derived from the *Obit class.
 *
 * This class has the facilities to convert from a GBT DCR file and
 * associated tables to OTF format.
 *
 * \section ObitGBTDCROTFUsage Usage
 * Instances can be obtained using the #newObitGBTDCROTF constructor,
 * the #ObitGBTDCROTFCopy constructor or a pointer duplicated using 
 * the #ObitGBTDCROTFRef macro.
 * When an instance is no longer needed, use the #ObitGBTDCROTFUnref 
 * macro to release it.
 *
 * To use this facility, Obit must have been started up defining a FITS 
 * directory and the base of the GBT archive structure, i.e. the directory
 * Antenna, DCR, and IF directories
 */

/*---------------Class Structure---------------------------*/
/** ObitGBTDCROTF Class. */
typedef struct {
#include "ObitGBTDCROTFDef.h"   /* actual definition */
} ObitGBTDCROTF;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitGBTDCROTF
 * returns a ObitGBTDCROTF*.
 * in = object to unreference
 */
#define ObitGBTDCROTFUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitGBTDCROTF.
 * returns a ObitGBTDCROTF*.
 * in = object to reference
 */
#define ObitGBTDCROTFRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitGBTDCROTFIsA(in) ObitIsA (in, ObitGBTDCROTFGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
 void ObitGBTDCROTFClassInit (void);

/** Public: Constructor. */
ObitGBTDCROTF* newObitGBTDCROTF (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitGBTDCROTFGetClass (void);

/** Public: Constructor from values. */
ObitGBTDCROTF* 
newObitGBTDCROTFValue (gchar *name, ObitOTF *outOTF, ObitErr *err);

/** Public: Convert GBT DCR files for a scan, append to OTF format */
ObitIOCode ObitGBTDCROTFConvert (ObitGBTDCROTF *in, olong inDisk, gchar *scanName, 
				 ObitErr *err);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * and function pointers.
 */
typedef struct  {
#include "ObitGBTDCROTFClassDef.h" /* Actual definition */
} ObitGBTDCROTFClassInfo; 


#endif /* OBITGBTDCR_H */ 
