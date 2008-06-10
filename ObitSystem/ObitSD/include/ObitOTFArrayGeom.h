/* $Id: ObitOTFArrayGeom.h,v 1.4 2005/10/06 19:33:28 bcotton Exp $ */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITOTFARRAYGEOM_H 
#define OBITOTFARRAYGEOM_H 

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glib.h>
#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitOTFSkyModel.h"
#include "ObitOTFDesc.h"
#include "ObitTableOTFArrayGeom.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFArrayGeom.h
 * GBT/OTF feed geometry
 *
 * This class is for creating and manipulating Array geometry description objects
 * for OTF data for the GBT
 * 
 * \section ObitOTFArrayGeomaccess Creators and Destructors
 * An ObitOTFArrayGeom will usually be created using ObitOTFArrayGeomCreate which allows 
 * specifying a name for the object as well as dimensionality of the array.
 *
 * A copy of a pointer to an ObitOTFArrayGeom should always be made using the
 * #ObitOTFArrayGeomRef function which updates the reference count in the object.
 * Then whenever freeing an ObitOTFArrayGeom or changing a pointer, the function
 * #ObitOTFArrayGeomUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitOTFArrayGeom Class structure. */
typedef struct {
#include "ObitOTFArrayGeomDef.h"   /* this class definition */
} ObitOTFArrayGeom;

/*----------------------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitOTFArrayGeom
 * returns a ObitOTFArrayGeom*.
 * in = object to unreference
 */
#define ObitOTFArrayGeomUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitOTFArrayGeom.
 * returns a ObitOTFArrayGeom*.
 * in = object to reference
 */
#define ObitOTFArrayGeomRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitOTFArrayGeomIsA(in) ObitIsA (in, ObitOTFArrayGeomGetClass())

/*---------------Public functions---------------------------*/
/**  Publi: Class initializer. */
void ObitOTFArrayGeomClassInit (void);

/** Public: Default Constructor. */
ObitOTFArrayGeom* newObitOTFArrayGeom (gchar* name);

/** Public: Create/initialize ObitOTFArrayGeom structures */
ObitOTFArrayGeom* ObitOTFArrayGeomCreate (olong ndetect);
/** Typedef for definition of class pointer structure */
typedef void (*ObitOTFArrayGeomCreateFP) (olong ndetect);

/** Public: ClassInfo pointer */
gconstpointer ObitOTFArrayGeomGetClass (void);

/** Public: Copy (deep) constructor. */
ObitOTFArrayGeom* 
ObitOTFArrayGeomCopy  (ObitOTFArrayGeom *in, ObitOTFArrayGeom *out, ObitErr *err);

/** Public: Copy (deep) constructor averaging over frequency . */
ObitOTFArrayGeom* 
ObitOTFArrayGeomAver  (ObitOTFArrayGeom *in, ObitOTFDesc *inDesc, 
		       ObitOTFArrayGeom *out, ObitOTFDesc *outDesc, 
		       ObitErr *err);

/** Public:  Read Table from disk */
ObitIOCode ObitOTFArrayGeomRead (ObitOTFArrayGeom **in, ObitTableOTFArrayGeom *table, ObitErr *err);

/** Public:  Write Table to disk */
ObitIOCode ObitOTFArrayGeomWrite (ObitOTFArrayGeom *in, ObitTableOTFArrayGeom *table, ObitErr *err);

/** Public: Get Parallactic angle for a given time and direction */
ofloat ObitOTFArrayGeomParAng (ObitOTFArrayGeom *in, ofloat time, ofloat ra, ofloat dec);

/** Public: Get Elevation for a given time and direction */
ofloat ObitOTFArrayGeomElev (ObitOTFArrayGeom *in, ofloat time, ofloat ra, ofloat dec);

/** Public: Get detector coordinates on sky */
void ObitOTFArrayGeomCoord(ObitOTFArrayGeom *in, ofloat raPoint, ofloat decPoint, ofloat rot,
			   ofloat *x, ofloat *y);

/** Public: Get detector locations projected onto a plane */
void ObitOTFArrayGeomProj(ObitOTFArrayGeom *in, ofloat raPoint, ofloat decPoint, ofloat rot,
			  ofloat raProj, ofloat decProj, ObitOTFProj Proj, ofloat *x, ofloat *y);

/** Public: Offset a celestial position in az, el */
void ObitOTFArrayGeomCorrPoint(ofloat azOff, ofloat elOff, ofloat pa,
			       ofloat *raPoint, ofloat *decPoint);
/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitOTFArrayGeomClassDef.h"
} ObitOTFArrayGeomClassInfo; 

#endif /* OBITOTFARRAYGEOM_H */ 
