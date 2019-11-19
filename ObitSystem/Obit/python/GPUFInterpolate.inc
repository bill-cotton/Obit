/* $Id:  $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for GPUFInterpolate type                   */
/*                                                                    */
/*;  Copyright (C) 2014,2019                                          */
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
#include "ObitGPUFInterpolate.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitGPUFInterpolate *me;
} GPUFInterpolate;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitGPUFInterpolate *me;
} GPUFInterpolate;

%inline %{
extern ObitGPUFInterpolate* 
GPUFInterpolateCreate (char* name, ObitFArray *inArray,  
                       ObitFArray *xArray, ObitFArray *yArray, 
		       long hwidth, ObitErr *err) {
  return ObitGPUFInterpolateCreate (name, inArray, xArray, yArray, hwidth, err);
} // end GPUFInterpolateCreate

ObitGPUFInterpolate* GPUFInterpolateRef (ObitGPUFInterpolate* in) {
  return ObitGPUFInterpolateRef (in);
} // end GPUFInterpolateRef

ObitGPUFInterpolate* GPUFInterpolateUnref (ObitGPUFInterpolate* in) {
  if (!ObitGPUFInterpolateIsA(in)) return NULL;
  return ObitGPUFInterpolateUnref (in);
} // end GPUFInterpolateUnref

void CreateGPUFInterpolate(PyObject *self, char* name, ObitFArray *inArray,  
                       ObitFArray *xArray, ObitFArray *yArray, 
		       long hwidth, ObitErr *err) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_GPUFInterpolate, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateGPUFInterpolate: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((GPUFInterpolate*)ptr)->me = 
     GPUFInterpolateCreate(name, inArray, xArray, yArray, hwidth, err);
   else  ((GPUFInterpolate*)ptr)->me = NULL;
   return;
 } // end CreateGPUFInterpolate

void DeleteGPUFInterpolate(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_GPUFInterpolate, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteGPUFInterpolate: could not recover c struct");
	return;
   }
   ((GPUFInterpolate*)ptr)->me = ObitGPUFInterpolateUnref(((GPUFInterpolate*)ptr)->me);
}// end DeleteGPUFInterpolate

ObitGPUFInterpolate* GPUFInterpolate_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_GPUFInterpolate, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"GPUFInterpolate_Get_me: could not recover c struct");
	return NULL;
   }
   return ((GPUFInterpolate*)ptr)->me;
}  // end GPUFInterpolate_Get_me

void GPUFInterpolate_Set_me(PyObject *self, ObitGPUFInterpolate *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_GPUFInterpolate, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"GPUFInterpolate_Set_me: could not recover c struct");
	return;
   }
   ((GPUFInterpolate*)ptr)->me =  ObitGPUFInterpolateRef(me);
}  // end GPUFInterpolate_Set_me

extern void 
GPUFInterpolateImage (ObitGPUFInterpolate *in, ObitFArray *inArray,
		      ObitFArray *outArray, ObitErr *err) {
  ObitGPUFInterpolateImage (in, inArray, outArray, err);
} // end GPUFInterpolateImage

extern ObitGPUFInterpolate* 
GPUFInterpolateCopy  (ObitGPUFInterpolate *in, ObitGPUFInterpolate *out, 
		      ObitErr *err) {
  return ObitGPUFInterpolateCopy (in, out, err);
} // end GPUFInterpolateCopy

extern void 
GPUFInterpolateClone (ObitGPUFInterpolate *in, ObitGPUFInterpolate *out, ObitErr *err) {
  ObitGPUFInterpolateClone (in, out, err);
} // end GPUFInterpolateClone 

extern long GPUFInterpolateIsA (ObitGPUFInterpolate* in) {
  return (long)ObitGPUFInterpolateIsA(in);
}

%}

