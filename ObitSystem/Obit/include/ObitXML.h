/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2008                                          */
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
#ifndef OBITXML_H 
#define OBITXML_H 

#include <xmlrpc.h>
#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitDConCleanWindow.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitXML.h
 * ObitXML XML class
 *
 * This class is derived from the #Obit class.
 *
 * This class controls the conversion to and from xml for 
 * Obit related information.  The xml structures are generally to be 
 * used in a remote procedure call using ObitRPC.
 * 
 * \section ObitXMLaccess Creators and Destructors
 * An ObitXML will usually be created using ObitXMLCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitXML should always be made using the
 * #ObitXMLRef function which updates the reference count in the object.
 * Then whenever freeing an ObitXML or changing a pointer, the function
 * #ObitXMLUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitXMLType
 * enum for object type.
 * This specifies which known type of XML object
 */
enum obitXMLType {
  /** Reply status from RPC call */
  OBIT_XML_Reply = 0, 
  /** Reply results */
  OBIT_XML_Result, 
  /** ObitInfoList */
  OBIT_XML_InfoList, 
  /** Ping parameters */
  OBIT_XML_Ping, 
  /** LoadImage parameters */
  OBIT_XML_LoadImage,
  /** EditWindow parameters */
  OBIT_XML_EditWindow,
  /** Binary blob */
  OBIT_XML_BinBlob,
  /** Unknown */
  OBIT_XML_Unknown
}; /* end enum obitXMLType */
/** typedef for enum for ObitXMLType object status. */
typedef enum obitXMLType ObitXMLType;


/*--------------Class definitions-------------------------------------*/
/** ObitXML Class structure. */
typedef struct {
#include "ObitXMLDef.h"   /* this class definition */
} ObitXML;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitXML
 * returns a ObitXML*.
 * in = object to unreference
 */
#define ObitXMLUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitXML.
 * returns a ObitXML*.
 * in = object to reference
 */
#define ObitXMLRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitXMLIsA(in) ObitIsA (in, ObitXMLGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitXMLClassInit (void);

/** Public: Default Constructor. */
ObitXML* newObitXML (gchar* name);

/** Public: ClassInfo pointer */
gconstpointer ObitXMLGetClass (void);

/** Public: Constructor for ping. */
ObitXML* 
ObitXMLPing2XML (ObitErr *err);

/** Public: Convert to ping info. */
gint
ObitXMLXML2Ping (ObitXML *xml, ObitErr *err);

/** Public: Constructor from ObitInfoList. */
ObitXML* 
ObitXMLInfoList2XML (ObitInfoList *list, ObitErr *err);

/** Public: Convert to ObitInfoList. */
ObitInfoList*
ObitXMLXML2InfoList (ObitXML *xml, ObitErr *err);

/** Public: Convert  XML created from an ObitInfoList to ObitInfoList. */
ObitInfoList*
ObitXMLXMLInfo2List (ObitXML *xml, ObitErr *err);

/** Public: Constructor from Image file info. */
ObitXML* 
ObitXMLFileInfo2XML (ObitIOType Type, gchar *Name,
		     gchar *AClass, gchar *ADir, olong ASeq, olong AUser,
		     olong Field, olong NField, 
		     ObitErr *err);

/** Public: Convert to Image file info. */
void
ObitXMLXML2FileInfo (ObitXML *xml, ObitIOType *Type, gchar **Name,
		     gchar **AClass, gchar **ADir, olong *ASeq, olong *AUser,
		     olong *Field, olong *NField, 
		     ObitErr *err);

/** Public: Constructor from ObitDConCleanWindow. */
ObitXML* 
ObitXMLWindow2XML (ObitDConCleanWindow *window, olong field, 
		   ObitErr *err);

/** Public: Convert to ObitDConCleanWindow. */
ObitDConCleanWindow* 
ObitXMLXML2Window (ObitXML *xml, ObitErr *err);

/** Public: Constructor from Binary blob */
ObitXML* 
ObitXMLBlob2XML (gpointer blob, ObitInfoList *desc, ObitErr *err);

/** Public: Convert to Binary blob. */
gpointer 
ObitXMLXML2Blob (ObitXML *xml, ObitInfoList **desc, ObitErr *err);

/** Public: Constructor for return object */
ObitXML* 
ObitXMLReturn (gchar *name, gpointer parmP,  ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitXMLClassDef.h"
} ObitXMLClassInfo; 

#endif /* OBITFXML_H */ 
