/* $Id$ */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITAIPSOBJECT_H 
#define OBITAIPSOBJECT_H 
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"
#include "ObitAIPS.h"
#include "ObitAIPSCat.h"

/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitAIPSObject.h
 * ObitAIPSObject module definition.
 *
 * This is a Utility module to handle the interface with the AIPS object
 * system.
 */

/** Maximum number of dimensions */
#define MAXAIPSOBJDIM 5

/*-------------- type definitions----------------------------------*/
/** Typedef for AIPS Object name */
typedef gchar AIPSObj[32];

/** Typedef for AIPS Object class */
typedef gchar AIPSObjClass[8];

/** Typedef for AIPS Object keyword */
typedef gchar AIPSKey[8];

/** Typedef for AIPS Object keyword dimension */
typedef oint AIPSKeyDim[MAXAIPSOBJDIM];

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitAIPSObjectType
 * enum for AIPS object keyword type. Defined in $INC/PAOOF.INC
 */
enum obitAIPSObjectType {
  /** Integer type */
  OBIT_AIPSObjectInt = 4,
  /** Float type */
  OBIT_AIPSObjectRe  = 2,
  /** Double type */
  OBIT_AIPSObjectDP  = 1,
   /** Character type */
  OBIT_AIPSObjectCar = 3,
   /** Boolean type */
  OBIT_AIPSObjectLog = 5,
}; /* end enum obitAIPSObjectType */

/** typedef for enum for AIPS Object keyword type */
typedef enum obitAIPSObjectType ObitAIPSObjectType;

/*-------------------- unions -------------------------------------*/
/** AIPS catalog header */
  union ObitAIPSCatEquiv { 
    oint   itg[256];
    ofloat flt[256];
    double dbl[256];
  };

/*---------------Public functions---------------------------*/
/** Public: Initialize AIPS Object manager. */
void ObitAIPSObjectOBinit (oint *ierr);

/** Public: Add a Catalog header keyword to the virtual keyword list */
void ObitAIPSObjectOBvhkw (AIPSObjClass class, AIPSKey keyword, ObitAIPSObjectType type, oint *ierr);

/** Public: See if keyword is an object dependent, virtual keyword. */
void ObitAIPSObjectOBkeyv (oint objnum, AIPSKey keywrd, oint *keypnt, olong *ierr, ObitErr *err);

/** Public: Fetch the value (array) for a specified real (non-virtual) keyword.*/
  void ObitAIPSObjectOBrget (oint objnum, AIPSKey keywrd, ObitAIPSObjectType *type, 
			     AIPSKeyDim dim, gpointer value, gchar *valuec, 
			     ObitErr *err);

/** Public: Associate an object slot with an object name. */
void ObitAIPSObjectOBcrea (AIPSObj name, AIPSObjClass class, ObitErr *err);

/** Public: Free the object slot associated with an object. */
void ObitAIPSObjectOBfree (AIPSObj name, ObitErr *err);

/** Public: Look up the object slot number of object with name "name". */
void  ObitAIPSObjectOBname (AIPSObj name, oint *objnum, ObitErr *err);

/** Public: Look up the class number/name of object number objnum. */
void  ObitAIPSObjectOBclass (oint objnum, oint *classno, AIPSObjClass name, ObitErr *err);

/** Public: Save an entry in an object creating it if necessary. */
void ObitAIPSObjectOBput (oint objnum, AIPSKey keywrd, ObitAIPSObjectType type, 
			  AIPSKeyDim dim, gpointer value, gchar *valuec, 
			  ObitErr *err);

/** Public: Fetch the value (array) for a specified keyword. */
void ObitAIPSObjectOBget (oint objnum, AIPSKey keywrd, ObitAIPSObjectType *type, 
			  AIPSKeyDim dim,  gpointer value, gchar *valuec, 
			  ObitErr *err);

/** Public: Check for for a specified keyword. */
gboolean 
ObitAIPSObjectOBinfo (oint objnum, AIPSKey keywrd, ObitAIPSObjectType *type, 
		      AIPSKeyDim dim, ObitErr *err);

/** Public: Return Disk and slot information for object. */
void ObitAIPSObjectOBdskc (AIPSObj name, oint *disk, oint *cno, ObitErr *err);

/** Public: Return catalog header record for an object. */
void ObitAIPSObjectOBhget (AIPSObj name, union ObitAIPSCatEquiv cat, ObitErr *err);

/** Public: Store catalog header record for an object. */
void ObitAIPSObjectOBhput (AIPSObj name, union ObitAIPSCatEquiv cat, ObitErr *err);

/** Public: Copies one image to another */
void ObitAIPSObjectOBcopy (AIPSObj namein, AIPSObj namout, ObitErr *err);

#endif /* OBITAIPSOBJECT_H */ 

