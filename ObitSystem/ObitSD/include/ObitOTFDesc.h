/* $Id: ObitOTFDesc.h,v 1.4 2006/03/22 21:23:20 bcotton Exp $    */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
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
#ifndef OBITOTFDESC_H 
#define OBITOTFDESC_H 
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"
#include "ObitImageDesc.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFDesc.h
 * ObitOTFDesc Obit "On the Fly" data descriptor class definition.
 * This class is derived from the Obit class.
 * This contains information about the observations and the size and 
 * structure of the data.
 *
 * \section ObitOTFDescUsage Usage
 * Instances can be obtained using the #newObitOTFDesc constructor
 * the #ObitOTFDescCopy copy constructor or a pointer duplicated using 
 * the #ObitOTFDescRef function.
 * When an instance is no longer needed, use the #ObitOTFDescUnref macro
 * to release it.
 */

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitOTFDesc
 * returns a ObitOTFDesc* (NULL).
 * \li in = object to unreference.
 */
#define ObitOTFDescUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitOTFDesc.
 * returns a ObitOTFDesc*.
 * in = object to reference
 */
#define ObitOTFDescRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitOTFDescIsA(in) ObitIsA (in, ObitOTFDescGetClass())

/** maximum data array dimension */
#define OTF_MAXDIM 7       
/** Maximum number of columns in data */
#define OTF_MAX_COL 10
/** Maximum length of descriptor string value */
#define OTFLEN_VALUE 41
/** Maximum length of descriptor keyword  */
#define OTFLEN_KEYWORD 21

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitGBTOTFType
 * enum for GBT OTF type
 * This specifies the GBT OTF data source 
 */
enum obitGBTOTFType {
  /** Unspecified  */
  OBIT_GBTOTF_Unknown = 0, 
  /** DCR */
  OBIT_GBTOTF_DCR,  
  /** Spectral processor */
  OBIT_GBTOTF_SP,
  /** CalTech Continuum Backend */
  OBIT_GBTOTF_CCB,
  /** Penn Array Receiver */
  OBIT_GBTOTF_PAR
}; /* end enum obitGBTOTFType */
/** typedef for enum for GBT OTF data source */
typedef enum obitGBTOTFType ObitGBTOTFType;

/*--------------Class definitions-------------------------------------*/
/**
 * ObitOTFDesc Class structure.
 *
 * This class contains descriptions of interferometric visibility data.
 */  
typedef struct {
#include "ObitOTFDescDef.h"  /* Actual definitions */
} ObitOTFDesc;

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitOTFDescClassInit (void);

/** Public: Constructor. */
ObitOTFDesc* newObitOTFDesc (gchar *name);

/** Public: Copy OTFDesc */
ObitOTFDesc* ObitOTFDescCopy (ObitOTFDesc* in, ObitOTFDesc* out,
			    ObitErr *err);

/** Public: Return class pointer. */
gconstpointer ObitOTFDescGetClass (void);

/** Public: Copy descriptive (nonstructural) information. */
void ObitOTFDescCopyDesc (ObitOTFDesc* in, ObitOTFDesc* out,
			 ObitErr *err);

/** Public: Index for easier access */
void ObitOTFDescIndex (ObitOTFDesc* in);

/** Public: Convert date string to JD */
void ObitOTFDescDate2JD (const gchar* date, odouble *JD);

/** Public: Convert JD to date string */
void ObitOTFDescJD2Date (odouble JD, gchar *date);

/** Public: OTF type to string */
void ObitOTFDescType2String (ObitGBTOTFType OTFType, gchar *TString);

/** Public: OTF string to type  */
ObitGBTOTFType ObitOTFDescString2Type (gchar *TString);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitOTFDescClassDef.h" /* Actual definition */
} ObitOTFDescClassInfo; 


#endif /* OBITOTFDESC_H */ 

