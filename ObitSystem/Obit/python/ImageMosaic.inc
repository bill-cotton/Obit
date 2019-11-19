/* $Id$   */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageMosaic type                       */
/*                                                                    */
/*;  Copyright (C) 2004-2019                                          */
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
#include "ObitImageMosaic.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitImageMosaic *me;
} ImageMosaic;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitImageMosaic *me;
} ImageMosaic;

%inline %{
extern ObitImageMosaic* newImageMosaic (char* name, long number) {
  return newObitImageMosaic (name, (olong)number);
} // end  newImageMosaic

extern ObitImageMosaic* ImageMosaicUnref (ObitImageMosaic* in) {
  if (!ObitImageMosaicIsA(in)) return NULL;
  return ObitImageMosaicUnref(in);
}

extern ObitImageMosaic*  ImageMosaicRef (ObitImageMosaic* in) {
  return ObitImageMosaicRef(in);
}

void CreateImageMosaic(PyObject *self, char* name, long number) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ImageMosaic, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateImageMosaic: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((ImageMosaic*)ptr)->me = newImageMosaic(name, number);
   else  ((ImageMosaic*)ptr)->me = NULL;
   return;
 } // end CreateImageMosaic

void DeleteImageMosaic(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ImageMosaic, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteImageMosaic: could not recover c struct");
	return;
   }
   ((ImageMosaic*)ptr)->me = ObitImageMosaicUnref(((ImageMosaic*)ptr)->me);
}// end DeleteImageMosaic

ObitImageMosaic* ImageMosaic_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ImageMosaic, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"ImageMosaic_Get_me: could not recover c struct");
	return NULL;
   }
   return ((ImageMosaic*)ptr)->me;
}  // end ImageMosaic_Get_me

void ImageMosaic_Set_me(PyObject *self, ObitImageMosaic *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ImageMosaic, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"ImageMosaic_Set_me: could not recover c struct");
	return;
   }
   ((ImageMosaic*)ptr)->me =  ObitImageMosaicRef(me);
}  // end ImageMosaic_Set_me

extern ObitImageMosaic* ImageMosaicCopy  (ObitImageMosaic *in, ObitImageMosaic *out, 
				          ObitErr *err) {
  return ObitImageMosaicCopy (in, out, err);
} // end  ImageMosaicCopy

extern void  ImageMosaicZapImage(ObitImageMosaic* in, long number, 
                                 ObitErr *err ) {
  return ObitImageMosaicZapImage(in, (olong)number, err);
}

extern ObitInfoList* ImageMosaicGetList (ObitImageMosaic* in) {
  return ObitInfoListRef(in->info);
}

extern ObitImage* ImageMosaicGetFullImage (ObitImageMosaic* in, ObitErr *err) {
  return ObitImageMosaicGetFullImage(in, err);
}

extern void ImageMosaicSetFullImage (ObitImageMosaic* in, 
                                 ObitImage *image, ObitErr *err) {
 ObitImageMosaicSetFullImage(in, image, err);
}

extern ObitImage* ImageMosaicGetImage (ObitImageMosaic* in, long number, 
                                       ObitErr *err) {
  return ObitImageMosaicGetImage(in, (olong)number, err);
}

extern void ImageMosaicSetImage (ObitImageMosaic* in, long number, 
                                 ObitImage *image, ObitErr *err) {
 ObitImageMosaicSetImage(in, (olong)number, image, err);
}

extern ObitImageMosaic* ImageMosaicCreate (char *name, ObitUV *uvData, ObitErr *err) {
 return ObitImageMosaicCreate(name, uvData, err);
}

extern void ImageMosaicDefine (ObitImageMosaic* in, ObitUV *uvData, long doBeam,
                               ObitErr *err) {
 gboolean LdoBeam = doBeam;
 ObitImageMosaicDefine(in, uvData, LdoBeam, err);
}

extern void ImageMosaicFlatten (ObitImageMosaic* in, ObitErr *err) {
 ObitImageMosaicFlatten(in, err);
}

extern char* ImageMosaicGetName (ObitImageMosaic* in) {
  if (ObitImageMosaicIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern long ImageMosaicGetNumber (ObitImageMosaic* in) {
  if (ObitImageMosaicIsA(in)) {
    return (long)in->numberImages;
  } else {
    return 0;
  }
}

extern long ImageMosaicIsA (ObitImageMosaic* in) {
  return (long)ObitImageMosaicIsA(in);
}

%}


