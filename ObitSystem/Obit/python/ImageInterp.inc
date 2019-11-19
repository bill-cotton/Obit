/* $Id: $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageInterp type                       */
/*                                                                    */
/*;  Copyright (C) 2009,2019                                          */
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
#include "ObitImageInterp.h"
#include "ObitImage.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitImageInterp *me;
} ImageInterp;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitImageInterp *me;
} ImageInterp;

%inline %{
extern ObitImageInterp* newImageInterp (char* name) {
  return newObitImageInterp (name);
} // end  newImageInterp

extern ObitImageInterp* ImageInterpUnref (ObitImageInterp* in) {
  if (!ObitImageInterpIsA(in)) return NULL;
  return ObitImageInterpUnref(in);
}

extern ObitImageInterp*  ImageInterpRef (ObitImageInterp* in) {
  return ObitImageInterpRef(in);
}

extern ObitImageInterp* ImageInterpCreate (char *name, ObitImage *image, 
				     long hwidth, ObitErr *err) {
  return ObitImageInterpCreate((gchar*)name, image, (olong)hwidth, err);
} // end ImageInterpCreate

void CreateImageInterp(PyObject *self, char* name, ObitImage *image, 
                       long hwidth, ObitErr *err) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ImageInterp, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateImageInterp: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((ImageInterp*)ptr)->me = ImageInterpCreate(name, image, hwidth, err);
   else  ((ImageInterp*)ptr)->me = NULL;
   return;
 } // end CreateImageInterp

void DeleteImageInterp(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ImageInterp, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteImageInterp: could not recover c struct");
	return;
   }
   ((ImageInterp*)ptr)->me = ObitImageInterpUnref(((ImageInterp*)ptr)->me);
}// end DeleteImageInterp

ObitImageInterp* ImageInterp_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ImageInterp, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"ImageInterp_Get_me: could not recover c struct");
	return NULL;
   }
   return ((ImageInterp*)ptr)->me;
}  // end ImageInterp_Get_me

void ImageInterp_Set_me(PyObject *self, ObitImageInterp *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ImageInterp, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"ImageInterp_Set_me: could not recover c struct");
	return;
   }
   ((ImageInterp*)ptr)->me =  ObitImageInterpRef(me);
}  // end ImageInterp_Set_me

extern ObitImageInterp* ImageInterpCopy  (ObitImageInterp *in, ObitImageInterp *out, 
				    ObitErr *err) {
  return ObitImageInterpCopy (in, out, err);
} // end  ImageInterpCopy

extern float ImageInterpValue (ObitImageInterp* in, double ra, double dec,
		            float parAng, long plane, ObitErr *err) {
  return ObitImageInterpValue(in, (odouble)ra, (odouble)dec, (ofloat)parAng, (olong)plane, err);
} // end ImageInterpValue

extern long ImageInterpFindPlane (ObitImageInterp* in, double freq) {
  return (long)ObitImageInterpFindPlane(in, (odouble)freq);
} // end ImageInterpFindPlane

extern char* ImageInterpGetName (ObitImageInterp* in) {
  if (ObitImageInterpIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern long ImageInterpIsA (ObitImageInterp* in) {
  return (long)ObitImageInterpIsA(in);
}

%}

