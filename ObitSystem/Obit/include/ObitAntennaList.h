/* $Id: ObitAntennaList.h,v 1.7 2007/09/03 15:34:30 bcotton Exp $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
#ifndef OBITANTENNALIST_H 
#define OBITANTENNALIST_H 
#include "Obit.h"
#include "ObitErr.h"
#include "ObitAntenna.h"
#include "ObitSource.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitAntennaList.h
 * ObitAntennaList class definition.
 *
 * This class is derived from the #Obit class.
 *
 * This class manages lists of antennas.
 *
 * \section ObitAntennaListUsage Usage
 * Instances can be obtained using the #newObitAntennaList constructor,
 * the #ObitAntennaListCopy constructor or a pointer duplicated using 
 * the #ObitAntennaListRef macro.
 * When an instance is no longer needed, use the #ObitAntennaListUnref 
 * macro to release it.
 */

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitIOStatus
 * enum for object status.
 * This specifies UV Polarization calibration type.
 */
enum obitUVPolCalType {
  /** No polarization calibration */
  OBIT_UVPoln_NoCal,
  /** Unrecognized code */
  OBIT_UVPoln_Unknown,
  /** R/L Linear D-term approximation */
  OBIT_UVPoln_Approx, 
  /** R/L Linear D-term approximation for resolved sources */
  OBIT_UVPoln_VLBI, 
  /** Elipticity-orientation */
  OBIT_UVPoln_ELORI,
  /** X/Y Linear D-term approximation */
  OBIT_UVPoln_XYLin
}; /* end enum obitIOStatus */
/** typedef for enum for ObitIO object status. */
typedef enum obitUVPolCalType ObitUVPolCalType;

/*---------------Class Structure---------------------------*/
/** ObitAntennaList Class. */
typedef struct {
#include "ObitAntennaListDef.h"   /* actual definition */
} ObitAntennaList;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitAntennaList
 * returns a ObitAntennaList*.
 * in = object to unreference
 */
#define ObitAntennaListUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitAntennaList.
 * returns a ObitAntennaList*.
 * in = object to reference
 */
#define ObitAntennaListRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitAntennaListIsA(in) ObitIsA (in, ObitAntennaListGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
 void ObitAntennaListClassInit (void);

/** Public: Constructor. */
ObitAntennaList* newObitAntennaList (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitAntennaListGetClass (void);

/** Public: Copy  constructor. */
ObitAntennaList* 
ObitAntennaListCopy  (ObitAntennaList *in, ObitAntennaList *out, ObitErr *err);

/** Public: Create from value */
ObitAntennaList* ObitAntennaListCreate (gchar* name, olong nant, olong numpolcal);

/** Public: Determine polarization calibration type */
ObitUVPolCalType ObitAntennaListGetPolType (gchar* type);

/** Public: Determine source elevation */
ofloat ObitAntennaListElev (ObitAntennaList *inAList, olong ant, ofloat time, 
			    ObitSource *Source);

/** Public: Determine source azimuth */
ofloat ObitAntennaListAz (ObitAntennaList *inAList, olong ant, ofloat time, 
			    ObitSource *Source);

/** Public: Determine antenna Parallactic angle */
ofloat ObitAntennaListParAng (ObitAntennaList *inAList, olong ant, ofloat time,
			      ObitSource *Source);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to parent class
 * and function pointers.
 */
typedef struct  {
#include "ObitAntennaListClassDef.h" /* Actual definition */
} ObitAntennaListClassInfo; 


#endif /* OBITANTENNALIST_H */ 
