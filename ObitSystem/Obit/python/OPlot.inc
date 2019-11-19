/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitPlot type                          */
/*                                                                    */
/*;  Copyright (C) 2006,2016,2019                                     */
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

%{
#include "ObitPlot.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitPlot *me;
} OPlot;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitPlot *me;
} OPlot;

%inline %{

/** Public: Create plot. */
extern ObitPlot* OPlotCreate(char* name) {
   return  newObitPlot ((gchar*)name);
}
extern ObitPlot* OPlotRef (ObitPlot* in) {
  return ObitPlotRef (in);
} // end OPlotRef

extern ObitPlot* OPlotUnref (ObitPlot* in) {
  if (!ObitPlotIsA(in)) return NULL;
  return ObitPlotUnref (in);
} // end OPlotUnref

void CreateOPlot(PyObject *self, char* name) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_OPlot, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Create: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((OPlot*)ptr)->me = OPlotCreate(name);
   else ((OPlot*)ptr)->me = NULL;
   return;
 } // end CreateOPlot

void DeleteOPlot(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_OPlot, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Delete: could not recover c struct");
	return;
   }
   ((OPlot*)ptr)->me = OPlotUnref(((OPlot*)ptr)->me);
} // end DeleteOPlot

ObitPlot* OPlot_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_OPlot, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Get_me: could not recover c struct");
	return NULL;
   }
   return ((OPlot*)ptr)->me;
} // end OPlot_Get_me

void OPlot_Set_me(PyObject *self, ObitPlot *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_OPlot, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"OPlot_Set_me: could not recover c struct");
	return;
   }
   ((OPlot*)ptr)->me =  OPlotRef(me);
}  // end OPlot_Set_me

/** Public: Initialize plot. */
extern void PlotInitPlot (ObitPlot* in, char *output, long color, 
                          long nx, long ny, ObitErr *err) {
  gchar *loutput=NULL;

  if (!strncmp (output, "None", 4)) {
    loutput = NULL;	
  } else {
    loutput = (gchar*)output;	
  }
  ObitPlotInitPlot (in, loutput, (olong)color, (olong)nx, (olong)ny, err);
} // end PlotInitPlot

/** Public: Finalize plot. */
extern void PlotFinishPlot (ObitPlot* in,  ObitErr *err) {
  ObitPlotFinishPlot (in, err);
} // end PlotFinishPlot

/** Public: Copy Plot */
extern ObitPlot* PlotCopy (ObitPlot* in, ObitPlot* out, ObitErr *err){
return  ObitPlotCopy (in, out, err);
} // end ObitPlotCopy

/** Public: Simple X-Y plot. */
extern void PlotXYPlot (ObitPlot* in, long symbol, long n, float *x, float *y, 
               	        ObitErr *err) {
  ObitPlotXYPlot (in, (olong)symbol, (olong)n, (ofloat*)x, (ofloat*)y, err);
} // end PlotXYPlot

/** Public: Simple X-Y over plot. */
extern void PlotXYOver (ObitPlot* in, long symbol, long n, float *x, float *y, 
	                ObitErr *err) {
  ObitPlotXYOver (in, (olong)symbol, (olong)n, (ofloat*)x, (ofloat*)y, err);
} // end PlotXYOver
		     
/** Public: Simple X-Y plot with error bars. */
extern void PlotXYErr (ObitPlot* in, long symbol, long n, float *x, float *y, 
	               float *e, ObitErr *err) {
  ObitPlotXYErr (in, (olong)symbol, (olong)n, (ofloat*)x, (ofloat*)y, (ofloat*)e, err);
} // end PlotXYErr


/** Public: Contour plot of image. */
extern void PlotContour (ObitPlot* in, char *label, ObitImage *image, float lev,
	                 float cntfac, ObitErr *err) {
  ObitPlotContour (in, (gchar*)label, image, (ofloat)lev, (ofloat)cntfac, err);
} // end PlotContour

/** Public: Gray Scale plot of image. */
extern void PlotGrayScale (ObitPlot* in, char *label, ObitImage *image, 
	                 ObitErr *err) {
  ObitPlotGrayScale (in, (gchar*)label, image, err);
} // end PlotGrayScale

/** Public: Mark positions on Contour plot of image. */
extern void PlotMarkCross (ObitPlot* in, ObitImage *image, long n,
		 	  double *ra, double *dec, float size, 
		          ObitErr *err) {
  ObitPlotMarkCross (in, image, (olong)n, (odouble*)ra, (odouble*)dec, (ofloat)size, err);
} // end PlotMarkCross


/** Public:  set window and viewport and draw labeled frame */
extern  void PlotSetPlot (ObitPlot* in, float xmin, float xmax, float ymin, float ymax, 
		         long just, long axis, ObitErr *err) {
   ObitPlotSetPlot (in, (ofloat)xmin, (ofloat)xmax, (ofloat)ymin, (ofloat)ymax, 
		    (olong)just, (olong)axis, err);
} // end PlotSetPlot

/** Public: write labels for x-axis, y-axis, and top of plot*/
extern  void PlotLabel (ObitPlot* in, char *xlabel, char *ylabel, char *title,
		        ObitErr *err)  {
  ObitPlotLabel (in, (gchar*)xlabel, (gchar*)ylabel, (gchar*)title,  err);
} // end ObitPlotLabel

/** Public: draw labeled frame around viewport */
extern void PlotDrawAxes (ObitPlot* in,  char *xopt, float xtick, long nxsub, 
		    char *yopt,  float ytick, long nysub, 
		    ObitErr *err) {
   ObitPlotDrawAxes (in, (gchar*)xopt, (ofloat)xtick, (olong)nxsub, 
		   (gchar*)yopt,  (ofloat)ytick, (olong)nysub, 
		    err);
} // end  PlotDrawAxes

/** Public: Scaling for characters */
extern  void PlotSetCharSize (ObitPlot* in, float cscale, ObitErr *err) {
   ObitPlotSetCharSize (in, (ofloat)cscale, err);
} // end PlotSetCharSize

/** Public: Set line width */
extern  void PlotSetLineWidth (ObitPlot* in, long lwidth, ObitErr *err) {
   ObitPlotSetLineWidth (in, (olong)lwidth, err);
} // end PlotSetLineWidth

/** Public: Set line style */
extern  void PlotSetLineStyle (ObitPlot* in, long lstyle, ObitErr *err) {
   ObitPlotSetLineStyle (in, (olong)lstyle, err);
} // end PlotSetLineStyle

/** Public: Set foreground color */
extern  void PlotSetColor (ObitPlot* in, long color, ObitErr *err) {
   ObitPlotSetColor (in, (olong)color, err);
} // end PlotSetColor

/** Public: Set/advance subpage */
extern  void PlotSetPage (ObitPlot* in, long sub, ObitErr *err) {
   ObitPlotSetPage (in, (olong)sub, err);
} // end PlotSetPage

/** Public: Write text */
extern void PlotText (ObitPlot* in, float x, float y,
		      float dx, float dy,  float fjust, char *text,
		      ObitErr *err) {
  ObitPlotText (in, (ofloat)x, (ofloat)y, (ofloat)dx, (ofloat)dy, 
		(ofloat)fjust, (gchar*)text,
	         err);
} // end PlotText

/** Public: Write text  relative to port */
extern void PlotRelText (ObitPlot* in, char *side, float disp, 
		      float coord, float fjust, char *text,
		       ObitErr *err) {
  ObitPlotRelText (in, (gchar*)side, (ofloat)disp, 
		(ofloat)coord, (ofloat)fjust, (gchar*)text,
	         err);
} // end PlotRelText

/**  Public: Draw a line..*/
extern void  PlotDrawLine (ObitPlot* in, float x1, float y1, 
                           float x2, float y2, ObitErr *err) {
   ObitPlotDrawLine (in, (ofloat)x1, (ofloat)y1, (ofloat)x2, (ofloat)y2, err);
} // end PlotDrawLine

/**  Public: Draw a curve.*/
extern void  PlotDrawCurve (ObitPlot* in, long n, float *x, float *y, 
                           ObitErr *err) {
   ObitPlotDrawCurve (in, (olong)n, (ofloat*)x, (ofloat*)y, err);
} // end PlotDrawCurve

/**  Public: Draw a circle.*/
extern void  PlotDrawCircle (ObitPlot* in, float x, float y, float radius, 
                           ObitErr *err) {
   ObitPlotDrawCircle (in, (ofloat)x, (ofloat)y, (ofloat)radius, err);
} // end PlotDrawCircle

/**  Public: Draw a Symbol.*/
extern void  PlotDrawSymbol (ObitPlot* in, float x, float y, long symbol, 
                             ObitErr *err) {
   ObitPlotDrawSymbol (in, (ofloat)x, (ofloat)y, (olong)symbol, err);
} // end PlotDrawSymbol

/**  Public: Draw a Polygon.*/
extern void  PlotDrawPoly (ObitPlot* in, long n, float *x, float *y, long fill, 
                           float scale, ObitErr *err) { 
   ObitPlotDrawPoly (in, (olong)n, (ofloat*)x, (ofloat*)y, (olong)fill, (ofloat)scale, err);
} // end PlotDrawPoly

extern long OPlotIsA (ObitPlot* in) {
  return (long)ObitPlotIsA(in);
} // end  OPlotIsA 

extern char* OPlotGetName (ObitPlot* in) {
  return in->name;
} // end  OPlotGetName

extern ObitInfoList* PlotGetList (ObitPlot* in) {
  return ObitInfoListRef(in->info);
}

%}

