/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageFit type                          */
/*                                                                    */
/*;  Copyright (C) 2007,2019                                          */
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
#include "ObitImageFit.h"
#include "ObitImage.h"
#include "ObitFitRegion.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitImageFit *me;
} ImageFit;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitImageFit *me;
} ImageFit;


%inline %{
extern ObitImageFit* newImageFit (char* name) {
  return newObitImageFit (name);
} // end  newImageFit

extern ObitImageFit* ImageFitUnref (ObitImageFit* in) {
  if (!ObitImageFitIsA(in)) return NULL;
  return ObitImageFitUnref(in);
}

extern ObitImageFit*  ImageFitRef (ObitImageFit* in) {
  return ObitImageFitRef(in);
}

void CreateImageFit(PyObject *self, char* name) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ImageFit, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateImageFit: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((ImageFit*)ptr)->me = newImageFit(name);
   else  ((ImageFit*)ptr)->me = NULL;
   return;
 } // end CreateImageFit

void DeleteImageFit(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ImageFit, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteImageFit: could not recover c struct");
	return;
   }
   ((ImageFit*)ptr)->me = ObitDConCleanUnref(((ImageFit*)ptr)->me);
} // end DeleteImageFit

ObitImageFit* ImageFit_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ImageFit, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"ImageFitGet_me: could not recover c struct");
	return NULL;
   }
   return ((ImageFit*)ptr)->me;
} // end ImageFit_Get_me

void ImageFit_Set_me(PyObject *self, ObitDConCleanWindow *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ImageFit, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"ImageFit_Set_me: could not recover c struct");
	return;
   }
   ((ImageFit*)ptr)->me =  ObitDConCleanWindowRef(me);
}  // end ImageFit_Set_me


extern ObitImageFit* ImageFitCopy  (ObitImageFit *in, ObitImageFit *out, 
				    ObitErr *err) {
  return ObitImageFitCopy (in, out, err);
} // end  ImageFitCopy

extern ObitInfoList* ImageFitGetList (ObitImageFit* in) {
  return ObitInfoListRef(in->info);
}

extern long ImageFitFit (ObitImageFit* in, ObitImage *image, 
                        ObitFitRegion* reg, ObitErr *err) {
 return (long)ObitImageFitFit(in, image, reg, err);
}

extern char* ImageFitGetName (ObitImageFit* in) {
  if (ObitImageFitIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern long ImageFitIsA (ObitImageFit* in) {
  return (long)ObitImageFitIsA(in);
}

%}

