/* $Id$    */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITPLOT_H 
#define OBITPLOT_H 
/* Disable plplot if not needed */
#ifdef HAVE_PLPLOT
#include <plplot.h>
/* Turn off HAVE_PGPLOT if plplot available */
#ifdef HAVE_PGPLOT
#undef HAVE_PGPLOT
#endif /* HAVE_ PGPLOT*/
#endif /* HAVE_PLPLOT*/
/* Disable pgplot if not needed */
#ifdef HAVE_PGPLOT
#include <cpgplot.h>
#endif /* HAVE_ PGPLOT*/
#include "Obit.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitImage.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPlot.h
 * ObitPlot Obit graphics class definition.
 *
 * This class is derived from the #Obit class.
 *
 * This contains information about the observations and the size and 
 * structure of the data.
 * Implementations using both pgplot and plplot with perference 
 * given to the second.
 *
 * \section ObitPlotUsage Usage
 * Instances can be obtained using the #newObitPlot constructor
 * the #ObitPlotCopy copy constructor or a pointer duplicated using 
 * the #ObitPlotRef function.
 * When an instance is no longer needed, use the #ObitPlotUnref macro
 * to release it.
 */

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitPlot
 * returns a ObitPlot* (NULL).
 * \li in = object to unreference.
 */
#define ObitPlotUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitPlot.
 * returns a ObitPlot*.
 * in = object to reference
 */
#define ObitPlotRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitPlotIsA(in) ObitIsA (in, ObitPlotGetClass())

/*--------------Class definitions-------------------------------------*/
/**
 * ObitPlot Class structure.
 *
 * This class contains descriptions of interferometric visibility data.
 */  
typedef struct {
#include "ObitPlotDef.h"  /* Actual definitions */
} ObitPlot;

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitPlotClassInit (void);

/** Public: Constructor. */
ObitPlot* newObitPlot (gchar *name);

/** Public: Initialize plot. */
void ObitPlotInitPlot (ObitPlot* in, gchar *output, olong color, 
		       olong nx, olong ny, ObitErr *err);

/** Public: Finalize plot. */
void ObitPlotFinishPlot (ObitPlot* in,  ObitErr *err);

/** Public: Copy Plot */
ObitPlot* ObitPlotCopy (ObitPlot* in, ObitPlot* out,
			    ObitErr *err);

/** Public: Return class pointer. */
gconstpointer ObitPlotGetClass (void);

/** Public: Simple X-Y plot. */
void ObitPlotXYPlot (ObitPlot* in, olong symbol, 
		     olong n, ofloat *x, ofloat *y, ObitErr *err);

/** Public: Simple X-Y over plot. */
void ObitPlotXYOver (ObitPlot* in, olong symbol, 
		     olong n, ofloat *x, ofloat *y, ObitErr *err);

/** Public: X-Y plot with error bars */
void ObitPlotXYErr (ObitPlot* in, olong symbol, 
		    olong n, ofloat *x, ofloat *y, ofloat *e, 
		    ObitErr *err);

/** Public: Contour plot */
void ObitPlotContour (ObitPlot* in, gchar *label, ObitImage *image,
		      ofloat lev, ofloat cntfac, ObitErr *err);

/** Public: Gray scale plot */
void ObitPlotGrayScale (ObitPlot* in, gchar *label, ObitImage *image,
			ObitErr *err);

/** Public: Mark positions on Contour plot */
void ObitPlotMarkCross (ObitPlot* in, ObitImage *image, olong n,
			odouble *ra, odouble *dec, ofloat size, ObitErr *err);

/** Public:  set window and viewport and draw labeled frame */
void  ObitPlotSetPlot (ObitPlot* in, 
		      ofloat xmin, ofloat xmax, ofloat ymin, ofloat ymax, 
		      olong just, olong axis, ObitErr *err);

/** Public: write labels for x-axis, y-axis, and top of plot*/
void  ObitPlotLabel (ObitPlot* in, gchar *xlabel, gchar *ylabel, gchar *title,
		      ObitErr *err);

/** Public: draw labeled frame around viewport */
void  ObitPlotDrawAxes (ObitPlot* in, 
		      gchar *xopt, ofloat xtick, olong nxsub, 
		      gchar *yopt,  ofloat ytick, olong nysub, 
		      ObitErr *err);

/** Public: Scaling for characters */
void  ObitPlotSetCharSize (ObitPlot* in, ofloat cscale, ObitErr *err);

/** Public: Set line width */
void  ObitPlotSetLineWidth (ObitPlot* in, olong lwidth, ObitErr *err);

/** Public: Set foreground color */
void  ObitPlotSetColor (ObitPlot* in, olong color, ObitErr *err);

/** Public: Set (sub)page */
void  ObitPlotSetPage (ObitPlot* in, olong sub, ObitErr *err);

/** Public: Write text */
void  ObitPlotText (ObitPlot* in, ofloat x, ofloat y, 
		    ofloat dx, ofloat dy,
		    ofloat just, gchar *text,
		    ObitErr *err);

/** Public: Write text relative to port */
void  ObitPlotRelText (ObitPlot* in, gchar *side, ofloat disp, 
		       ofloat coord, ofloat fjust, gchar *text,
		       ObitErr *err);

/**  Public: Draw a line.*/
void  ObitPlotDrawLine (ObitPlot* in,  ofloat x1, ofloat y1, 
			ofloat x2, ofloat y2, ObitErr *err);

/**  Public: Draw a curve.*/
void  ObitPlotDrawCurve (ObitPlot* in,  olong n, ofloat *x, ofloat *y, 
			 ObitErr *err);

/**  Public: Draw a Symbol.*/
void  ObitPlotDrawSymbol (ObitPlot* in,  ofloat x, ofloat y, 
			  olong symbol, ObitErr *err);

/*-------------------Class Info--------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitPlotClassDef.h" /* Actual definition */
} ObitPlotClassInfo; 


#endif /* OBITPLOT_H */ 

