/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2009                                          */
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
#ifndef OBITDCONCLEANWINDOW_H 
#define OBITDCONCLEANWINDOW_H 

#include "Obit.h"
#include "ObitImageMosaic.h"
#include "ObitErr.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanWindow.h
 * ObitDConCleanWindow class defining CLEAN windows.
 *
 * This class is derived from the #Obit class.
 *
 * This class contains specifications of which pixels in the images of an
 * #ObitImageMosaic are candidates for a CLEAN component.
 * The current implementation uses a GList for each field so that the number of windows
 * is arbitrary.  However, this should be transparent outside the class.
 * Each field in an ImageMosaic has two potential sets of CLEAN windows, 
 * the traditional "inner" window with an arbitrary number of components and an
 * "outer", single window which sets the region in which the autoWindow feature
 * is allowed to place windows.
 * "Unwindows" are inner windows in which CLEANing is NOT allowed.
 * 
 * \section ObitDConCleanWindowaccess Creators and Destructors
 * An ObitDConCleanWindow will usually be created using ObitDConCleanWindowCreate 
 * which allows specifying a name for the object as well as dimensionality of the array.
 *
 * A copy of a pointer to an ObitDConCleanWindow should always be made using the
 * #ObitDConCleanWindowRef function which updates the reference count in the object.
 * Then whenever freeing an ObitDConCleanWindow or changing a pointer, the function
 * #ObitDConCleanWindowUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 *
 * The contents can be modified using the #ObitDConCleanWindowAdd and 
 * #ObitDConCleanWindowUpdate functions.
 * The presence of valid pixels in an image is indicated by
 * #ObitDConCleanWindowImage and in a row by #ObitDConCleanWindowRow this also 
 * computes a mask for the row indicating valid pixels.
 *
 * This class supports the autoWindow facility and has different behavior when the 
 * autoWindow member is TRUE or FALSE.
 * If TRUE then if no inner window is specified then all pixels are invalid.
 * If FALSE then the default is that all pixels are selected.
 */

/*--------------Class definitions-------------------------------------*/
/** ObitDConCleanWindow Class structure. */
typedef struct {
#include "ObitDConCleanWindowDef.h"   /* this class definition */
} ObitDConCleanWindow;

/*----------------- enums ---------------------------*/
/**
 * enum for window types
 * Should be coordinated with ObitErrorLevelString in ObitErr.c.
 */
enum ObitDConCleanWindowType {
  /** rectangle, specified by blc, trc corners */
  OBIT_DConCleanWindow_rectangle = 0,
  /** round, specified by radius, center x, center y pixel */
  OBIT_DConCleanWindow_round,
  /** rectangle unwindow, specified by blc, trc corners */
  OBIT_DConCleanWindow_unrectangle,
  /** round unwindow, specified by radius, center x, center y pixel */
  OBIT_DConCleanWindow_unround
};/* end enum ObitDConObitCleanWindowType */
/**
 * typedef for enum for window types
 */
typedef enum ObitDConCleanWindowType ObitDConCleanWindowType;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitDConCleanWindow
 * returns a ObitDConCleanWindow*.
 * in = object to unreference
 */
#define ObitDConCleanWindowUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitDConCleanWindow.
 * returns a ObitDConCleanWindow*.
 * in = object to reference
 */
#define ObitDConCleanWindowRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitDConCleanWindowIsA(in) ObitIsA (in, ObitDConCleanWindowGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitDConCleanWindowClassInit (void);

/** Public: Default Constructor. */
ObitDConCleanWindow* newObitDConCleanWindow (gchar* name);

/** Public: Create/initialize ObitDConCleanWindow structures */
ObitDConCleanWindow* 
ObitDConCleanWindowCreate (gchar* name, ObitImageMosaic *mosaic, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitDConCleanWindow* (*ObitDConCleanWindowCreateFP) (gchar* name, 
							     ObitImageMosaic *mosaic, 
							     ObitErr *err);

/** Public: Create/initialize ObitDConCleanWindow structure with 1 field */
ObitDConCleanWindow* 
ObitDConCleanWindowCreate1 (gchar* name, olong naxis[2], ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitDConCleanWindow* 
(*ObitDConCleanWindowCreate1FP) (gchar* name, olong naxis[2], ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitDConCleanWindowGetClass (void);

/** Public: Copy (deep) constructor. */
ObitDConCleanWindow* ObitDConCleanWindowCopy  (ObitDConCleanWindow *in, 
					       ObitDConCleanWindow *out, 
					       ObitErr *err);

/** Public: Copy structure. */
void ObitDConCleanWindowClone (ObitDConCleanWindow *in, 
			       ObitDConCleanWindow *out, 
			       ObitErr *err);

/** Public: Ask window definition */
gboolean ObitDConCleanWindowInfo (ObitDConCleanWindow *in, 
			      olong field, olong Id,
			      ObitDConCleanWindowType *type,
			      olong **window,  ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef gboolean (*ObitDConCleanWindowInfoFP) (ObitDConCleanWindow *in, 
					       olong field, olong Id,
					       ObitDConCleanWindowType *type,
					       olong **window, ObitErr *err);

/** Public: Search for a window near a given pixel */
olong ObitDConCleanWindowSearch (ObitDConCleanWindow *in, 
				 olong field, olong pixel[2], 
				 olong toler, olong *which, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef olong (*ObitDConCleanWindowSearchFP) (ObitDConCleanWindow *in, 
					      olong field, olong pixel[2], 
					      olong toler, olong *which, 
					      ObitErr *err);

/** Public: Add a new window definition */
olong ObitDConCleanWindowAdd (ObitDConCleanWindow *in, 
			      olong field, ObitDConCleanWindowType type,
			      olong *window, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef olong (*ObitDConCleanWindowAddFP) (ObitDConCleanWindow *in, 
					   olong field, 
					   ObitDConCleanWindowType type,
					   olong *window, ObitErr *err);

/** Public: Delete a window */
void ObitDConCleanWindowDel (ObitDConCleanWindow *in, 
			     olong field, olong Id, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitDConCleanWindowDelFP) (ObitDConCleanWindow *in, 
					  olong field, olong Id, ObitErr *err);

/** Public: Modify an existing window */
void ObitDConCleanWindowUpdate (ObitDConCleanWindow *in,  
				olong field, olong Id, 
				ObitDConCleanWindowType type,
				olong *window, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitDConCleanWindowUpdateFP) (ObitDConCleanWindow *in,  
				olong field, olong Id, 
				ObitDConCleanWindowType type,
				olong *window, ObitErr *err);

/** Public: Set outer window for a field  */
void ObitDConCleanWindowOuter (ObitDConCleanWindow *in, olong field, 
			       ObitDConCleanWindowType type,
			       olong *window, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void (*ObitDConCleanWindowOuterFP) (ObitDConCleanWindow *in, 
					    olong field,  
					     ObitDConCleanWindowType type,
					     olong *window, ObitErr *err);

/** Public: Are there any valid pixels in this field's image? */
gboolean ObitDConCleanWindowImage (ObitDConCleanWindow *in, olong field, 
				   ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef gboolean (*ObitDConCleanWindowImageFP) (ObitDConCleanWindow *in, 
						olong field, ObitErr *err);

/** Public: Are there any valid pixels in a specified row within inner window?  */
gboolean ObitDConCleanWindowRow (ObitDConCleanWindow *in, olong field, olong row, 
				 gboolean **mask, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef gboolean (*ObitDConCleanWindowRowFP) (ObitDConCleanWindow *in, 
					      olong field, olong row, 
					      gboolean **mask, ObitErr *err);

/** Public: Are there any valid pixels in a specified row in positive boxes?  */
gboolean ObitDConCleanWindowInnerRow (ObitDConCleanWindow *in, olong field, olong row, 
				      gboolean **mask, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef gboolean (*ObitDConCleanWindowInnerRowFP) (ObitDConCleanWindow *in, 
						   olong field, olong row, 
						   gboolean **mask, ObitErr *err);

/** Public: Are there any valid pixels in a specified row in unboxes?  */
gboolean ObitDConCleanWindowUnrow (ObitDConCleanWindow *in, olong field, olong row, 
				 gboolean **mask, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef gboolean (*ObitDConCleanWindowUnrowFP) (ObitDConCleanWindow *in, 
					      olong field, olong row, 
					      gboolean **mask, ObitErr *err);

/** Public: Are there any valid pixels in a specified row with outer window?  */
gboolean ObitDConCleanWindowOuterRow (ObitDConCleanWindow *in, olong field, 
				      olong row, gboolean **mask, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef gboolean (*ObitDConCleanWindowOuterRowFP) (ObitDConCleanWindow *in, 
						   olong field, olong row, 
						   gboolean **mask, ObitErr *err);

/** Public: What is the maximum region covered in x or y?  */
olong ObitDConCleanWindowSize (ObitDConCleanWindow *in, olong field, 
				 ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef olong (*ObitDConCleanWindowSizeFP) (ObitDConCleanWindow *in, 
					      olong field, ObitErr *err);

/** Public: How many pixels are selected  */
olong ObitDConCleanWindowCount (ObitDConCleanWindow *in, olong field, 
				ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef olong (*ObitDConCleanWindowCountFP) (ObitDConCleanWindow *in, 
					     olong field, 
					     ObitDConCleanWindowType type,
					     olong *window, ObitErr *err);

/** Public: find values needed for autoWindow  */
gboolean 
ObitDConCleanWindowAutoWindow (ObitDConCleanWindow *in, 
			       olong field, ObitFArray *image,
			       gboolean doAbs,
			       ofloat *PeakIn, olong *PeakInPos,
			       ofloat *PeakOut, ofloat *RMS,
			       ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef gboolean 
(*ObitDConCleanWindowAutoWindowFP) (ObitDConCleanWindow *in, 
				    olong field, ObitFArray *image,
				    ofloat *PeakIn, olong *PeakInPos,
				    ofloat *PeakOut, ofloat *RMS,
				    ObitErr *err);

/** Public: Replace all windows for a given field with those from another window  */
void 
ObitDConCleanWindowReplaceField (ObitDConCleanWindow *in,  olong ifield, 
				 ObitDConCleanWindow *out, olong ofield,
				 ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef void 
(*ObitDConCleanWindowReplaceFieldFP) (ObitDConCleanWindow *in,  olong ifield, 
				 ObitDConCleanWindow *out, olong ofield,
				 ObitErr *err);

/** Public: Add a field to a window object */
olong 
ObitDConCleanWindowAddField (ObitDConCleanWindow *in,  
			     olong inaxes[2], ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef olong 
(*ObitDConCleanWindowAddFieldFP) (ObitDConCleanWindow *in,  
				  olong inaxes[2], ObitErr *err);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitDConCleanWindowClassDef.h"
} ObitDConCleanWindowClassInfo; 

#endif /* OBITFDCONCLEANWINDOW_H */ 
