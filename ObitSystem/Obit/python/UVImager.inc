/* $Id$ */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageMosaic type                       */
/*                                                                    */
/*;  Copyright (C) 2005,2009,2019                                     */
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
#include "ObitUVImager.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitUVImager *me;
} UVImager;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitUVImager *me;
} UVImager;


%inline %{
extern ObitUVImager* newUVImager (char* name) {
  return newObitUVImager (name);
} // end  newUVImager

extern ObitUVImager* UVImagerCreate (char *name, ObitUV *uvData, ObitErr *err) {
 return ObitUVImagerCreate(name, uvData, err);
}

extern ObitUVImager* UVImagerUnref (ObitUVImager* in) {
  if (!ObitUVImagerIsA(in)) return NULL;
  return ObitUVImagerUnref(in);
}

extern ObitUVImager*  UVImagerRef (ObitUVImager* in) {
  return ObitUVImagerRef(in);
}

void CreateUVImager(PyObject *self, char* name, ObitUV *uvData, ObitErr *err) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UVImager, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Create: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((UVImager*)ptr)->me = UVImagerCreate(name, uvData, err);
   else ((UVImager*)ptr)->me = NULL;
   return;
 } // end CreateUVImager

void DeleteUVImager(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UVImager, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Delete: could not recover c struct");
	return;
   }
   ((UVImager*)ptr)->me = UVImagerUnref(((UVImager*)ptr)->me);
} // end DeleteUVImager

ObitUVImager* UVImager_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UVImager, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Get_me: could not recover c struct");
	return NULL;
   }
   return ((UVImager*)ptr)->me;
} // end UVImager_Get_me

void UVImager_Set_me(PyObject *self, ObitUVImager *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UVImager, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"UVImager_Set_me: could not recover c struct");
	return;
   }
   ((UVImager*)ptr)->me =  UVImagerRef(me);
}  // end UVImager_Set_me

extern ObitUVImager* UVImagerCopy  (ObitUVImager *in, ObitUVImager *out, 
				    ObitErr *err) {
  return ObitUVImagerCopy (in, out, err);
} // end  UVImagerCopy

extern ObitUV* UVImagerGetUV (ObitUVImager* in) {
  return ObitUVRef(in->uvdata);
}

extern ObitImageMosaic* UVImagerGetMosaic (ObitUVImager* in, ObitErr *err) {
  return ObitUVImagerGetMosaic(in, err);
}

extern void UVImagerWeight (ObitUVImager* in, ObitErr *err) {
 ObitUVImagerWeight(in, err);
}

extern void UVImagerImage (ObitUVImager* in, long field, long doWeight, long doBeam, 
                            long doFlatten, ObitErr *err) {
 gboolean LdoWeight=doWeight, LdoBeam=doBeam, LdoFlatten=doFlatten;
 olong Lfield[2];
 Lfield[0] = field; Lfield[1] =  0;
 ObitUVImagerImage(in, Lfield, LdoWeight, LdoBeam, LdoFlatten, err);
}

extern void UVImagerFlatten (ObitUVImager* in, ObitErr *err) {
 ObitUVImagerFlatten(in, err);
}

extern char* UVImagerGetName (ObitUVImager* in) {
  if (ObitUVImagerIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern long UVImagerIsA (ObitUVImager* in) {
  return (long)ObitUVImagerIsA(in);
}

%}

