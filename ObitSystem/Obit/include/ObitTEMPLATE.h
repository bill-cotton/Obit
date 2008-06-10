/* $Id: ObitTEMPLATE.h,v 1.9 2007/08/31 17:24:48 bcotton Exp $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2007                                               */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/* CHANGE HERE TO REFLECT CLASS NAME                                  */
#ifndef OBITTEMPLATE_H 
#define OBITTEMPLATE_H 

/* CHANGE ?Def.h and ?ClassDef.h if parent class is other than Obit  */
#include "Obit.h"
#include "ObitErr.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/* CHANGE HERE                                                        */
/**
 * \file ObitTEMPLATE.h
 *
 * ObitTEMPLATE template for classes derived from #Obit
 *
 * Class documentation should go here.
 * 
 * \section ObitTEMPLATEaccess Creators and Destructors
 * An ObitTEMPLATE will usually be created using ObitTEMPLATECreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitTEMPLATE should always be made using the
 * #ObitTEMPLATERef function which updates the reference count in the object.
 * Then whenever freeing an ObitTEMPLATE or changing a pointer, the function
 * #ObitTEMPLATEUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitTEMPLATE Class structure. */
typedef struct {
#include "ObitTEMPLATEDef.h"   /* this class definition */
} ObitTEMPLATE;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitTEMPLATE
 * returns a ObitTEMPLATE*.
 * in = object to unreference
 */
#define ObitTEMPLATEUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitTEMPLATE.
 * returns a ObitTEMPLATE*.
 * in = object to reference
 */
#define ObitTEMPLATERef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitTEMPLATEIsA(in) ObitIsA (in, ObitTEMPLATEGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitTEMPLATEClassInit (void);

/** Public: Default Constructor. */
ObitTEMPLATE* newObitTEMPLATE (gchar* name);

/** Public: Create/initialize ObitTEMPLATE structures */
ObitTEMPLATE* ObitTEMPLATECreate (gchar* name);
/** Typedef for definition of class pointer structure */
typedef ObitTEMPLATE* (*ObitTEMPLATECreateFP) (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitTEMPLATEGetClass (void);

/** Public: Copy (deep) constructor. */
ObitTEMPLATE* ObitTEMPLATECopy  (ObitTEMPLATE *in, ObitTEMPLATE *out, ObitErr *err);

/** Public: Copy structure. */
void ObitTEMPLATEClone (ObitTEMPLATE *in, ObitTEMPLATE *out, ObitErr *err);

/* CHANGE HERE                                                        */


/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitTEMPLATEClassDef.h"
} ObitTEMPLATEClassInfo; 

#endif /* OBITFTEMPLATE_H */ 
