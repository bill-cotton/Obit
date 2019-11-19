/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for CleanVis type                          */
/*                                                                    */
/*;  Copyright (C) 2005-2019                                          */
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

%{
#include "ObitDConClean.h"
#include "ObitDConCleanVis.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitDConCleanVis *me;
} CleanVis;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitDConCleanVis *me;
} CleanVis;

%inline %{
extern ObitDConCleanVis* newCleanVis (char* name) {
  return newObitDConCleanVis (name);
} // end  newCleanVis

extern ObitDConCleanVis* CleanVisCreate (char *name, ObitUV* uvdata, 
                                         ObitErr *err) {
 return ObitDConCleanVisCreate(name, uvdata, err);
}

extern ObitDConCleanVis* CleanVisUnref (ObitDConCleanVis* in) {
  if (!ObitDConCleanVisIsA(in)) return NULL;
  return ObitDConCleanVisUnref(in);
}

extern ObitDConCleanVis*  CleanVisRef (ObitDConCleanVis* in) {
  return ObitDConCleanVisRef(in);
}

void CreateCleanVis(PyObject *self, char* name) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_CleanVis, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateCleanVis: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((CleanVis*)ptr)->me = newCleanVis(name);
   else  ((CleanVis*)ptr)->me = NULL;
   return;
 } // end CreateCleanVis

void DeleteCleanVis(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_CleanVis, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteCleanVis: could not recover c struct");
	return;
   }
   ((CleanVis*)ptr)->me = CleanVisUnref(((CleanVis*)ptr)->me);
}// end DeleteCleanVis

ObitDConCleanVis* CleanVis_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_CleanVis, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CleanVis_Get_me: could not recover c struct");
	return NULL;
   }
   return ((CleanVis*)ptr)->me;
}  // end CleanVis_Get_me

void CleanVis_Set_me(PyObject *self, ObitDConCleanVis *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_CleanVis, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CleanVis_Set_me: could not recover c struct");
	return;
   }
   ((CleanVis*)ptr)->me =  CleanVisRef(me);
}  // end CleanVis_Set_me

extern ObitDConCleanVis* CleanVisCopy  (ObitDConCleanVis *in, 
	ObitDConCleanVis *out, 
        ObitErr *err) {
  return ObitDConCleanVisCopy (in, out, err);
} // end  CleanVisCopy

extern ObitInfoList* CleanVisGetList (ObitDConCleanVis* in) {
  return ObitInfoListRef(in->info);
}

extern ObitImageMosaic* CleanVisGetImageMosaic (ObitDConCleanVis* in) {
  return ObitImageMosaicRef(in->mosaic);
}

extern void CleanVisSetImageMosaic (ObitDConCleanVis* in, ObitImageMosaic *mosaic) {
  in->mosaic = ObitImageMosaicUnref(in->mosaic);  /* Out with the old */
  in->mosaic = ObitImageMosaicRef(mosaic);        /* In with the new */
}

extern ObitSkyModel* CleanVisGetSkyModel (ObitDConCleanVis* in) {
  return ObitSkyModelRef(in->skyModel);
}

extern void CleanVisSetSkyModel (ObitDConCleanVis* in, ObitSkyModel *skymodel) {
  in->skyModel = ObitSkyModelUnref(in->skyModel);  /* Out with the old */
  in->skyModel = ObitSkyModelRef(skymodel);        /* In with the new */
}

extern ObitDConCleanWindow* CleanVisGetWindow (ObitDConCleanVis* in) {
  return ObitDConCleanWindowRef(in->window);
}

extern void CleanVisSetWindow (ObitDConCleanVis* in, ObitDConCleanWindow *window) {
  in->window = ObitDConCleanWindowUnref(in->window);  /* Out with the old */
  in->window = ObitDConCleanWindowRef(window);        /* In with the new */
}

// if (win[0]<0 then this is a round window then win[1]=radius, [2,3] = center
// else rectangular and blc=(win[0],win[1]), trc= blc=(win[2],win[3])
extern void CleanVisAddWindow (ObitDConCleanVis* in, long field, long *win, 
                                 ObitErr *err) {
  olong lfield, window[4];
  ObitDConCleanWindowType type;

  // which type? 
  if (win[0]<0) { // circle
    type = OBIT_DConCleanWindow_round;
    window[0] = win[1];
    window[1] = win[2];
    window[2] = win[3];
  } else { // rectangular
    type = OBIT_DConCleanWindow_rectangle;
    window[0] = win[0];
    window[1] = win[1];
    window[2] = win[2];
    window[3] = win[3];
  }
	
  lfield = field;
  ObitDConCleanWindowAdd (in->window, lfield, type, window, err);
}  // end CleanVisAddWindow 

extern void CleanVisDeconvolve (ObitDConCleanVis* in, ObitErr *err) {
 ObitDConCleanVisDeconvolve((ObitDCon*)in, err);
}

extern void CleanVisRestore (ObitDConCleanVis* in, ObitErr *err) {
 ObitDConCleanRestore ((ObitDConClean*)in, err);
}

extern void CleanVisDefWindow (ObitDConCleanVis* in, ObitErr *err) {
 ObitDConCleanVisDefWindow((ObitDConClean*)in, err);
}

extern long CleanVisReimage (ObitDConCleanVis* in, ObitUV* uvdata, 
                            ObitErr *err) {
  gboolean bout;
  long out;
  bout = ObitDConCleanVisReimage(in, uvdata, err);
  if (bout) out = 1;
  else out = 0;
  return out;
}

extern char* CleanVisGetName (ObitDConCleanVis* in) {
  if (ObitDConCleanVisIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern long CleanVisIsA (ObitDConCleanVis* in) {
  return (long)ObitDConCleanVisIsA(in);
}

%}

