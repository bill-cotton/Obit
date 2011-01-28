/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009                                               */
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
#ifndef OBITPRINTER_H 
#define OBITPRINTER_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitFile.h"
#include "ObitThread.h"
#include "ObitInfoList.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPrinter.h
 *
 * ObitPrinter template for classes derived from #Obit
 *
 * The printer class handles the actual printing, either interactively 
 * or to a file.
 * 
 * \section ObitPrinteraccess Creators and Destructors
 * An ObitPrinter will usually be created using ObitPrinterCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitPrinter should always be made using the
 * #ObitPrinterRef function which updates the reference count in the object.
 * Then whenever freeing an ObitPrinter or changing a pointer, the function
 * #ObitPrinterUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitPrinter Class structure. */
typedef struct {
#include "ObitPrinterDef.h"   /* this class definition */
} ObitPrinter;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitPrinter
 * returns a ObitPrinter*.
 * in = object to unreference
 */
#define ObitPrinterUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitPrinter.
 * returns a ObitPrinter*.
 * in = object to reference
 */
#define ObitPrinterRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitPrinterIsA(in) ObitIsA (in, ObitPrinterGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitPrinterClassInit (void);

/** Public: Default Constructor. */
ObitPrinter* newObitPrinter (gchar* name);

/** Public: Create/initialize ObitPrinter structures */
ObitPrinter* ObitPrinterCreate (gchar* name, gboolean isInteractive,
				FILE *outStream, gchar *fileName);
/** Typedef for definition of class pointer structure */
typedef ObitPrinter* (*ObitPrinterCreateFP) (gchar* name, 
					     gboolean isInteractive,
					     FILE *outStream, 
					     gchar *fileName);

/** Public: ClassInfo pointer */
gconstpointer ObitPrinterGetClass (void);

/** Public: Copy (deep) constructor. */
ObitPrinter* ObitPrinterCopy  (ObitPrinter *in, ObitPrinter *out, ObitErr *err);

/** Public: Copy structure. */
void ObitPrinterClone (ObitPrinter *in, ObitPrinter *out, ObitErr *err);

/** Public: Open Printer. */
void ObitPrinterOpen (ObitPrinter *in, olong LinesPerPage, gchar *Title1, 
		      gchar *Title2, ObitErr *err);
typedef void* (*ObitPrinterOpenFP) (ObitPrinter *in, olong LinesPerPage, 
				    gchar *Title1, gchar *Title2, ObitErr *err);

/** Public: Print line. */
void ObitPrinterWrite (ObitPrinter *in, gchar *line, gboolean *quit, 
		       ObitErr *err);
typedef void* (*ObitPrinterWriteFP) (ObitPrinter *in, gchar *line, 
				     gboolean *quit, ObitErr *err);

/** Public: Close Printer. */
void ObitPrinterClose (ObitPrinter *in, ObitErr *err);
typedef void* (*ObitPrinterCloseFP) (ObitPrinter *in, ObitErr *err);

/** Public: Update Titles. */
void ObitPrinterSetTitle (ObitPrinter *in, gchar *Title1, 
			  gchar *Title2, ObitErr *err);
typedef void* (*ObitPrinterSetTitleFP) (ObitPrinter *in, gchar *Title1, 
					gchar *Title2, ObitErr *err);

/** Public: Start a new page. */
void ObitPrinterNewPage (ObitPrinter *in, gboolean *quit, ObitErr *err);
typedef void* (*ObitPrinterNewPageFP) (ObitPrinter *in, gboolean *quit, 
				       ObitErr *err);


/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitPrinterClassDef.h"
} ObitPrinterClassInfo; 

#endif /* OBITFPRINTER_H */ 
