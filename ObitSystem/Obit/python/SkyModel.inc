/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for SkyModel type                          */
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
#include "ObitSkyModel.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitSkyModel *me;
} SkyModel;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitSkyModel *me;
} SkyModel;


%inline %{
extern ObitSkyModel* newSkyModel (char* name) {
  return newObitSkyModel (name);
} // end  newSkyModel

extern ObitSkyModel* SkyModelUnref (ObitSkyModel* in) {
  if (!ObitSkyModelIsA(in)) return NULL;
  return ObitSkyModelUnref(in);
}

extern ObitSkyModel*  SkyModelRef (ObitSkyModel* in) {
  return ObitSkyModelRef(in);
}

void CreateSkyModel(PyObject *self, char* name) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_SkyModel, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Create: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((SkyModel*)ptr)->me = newSkyModel(name);
   else ((SkyModel*)ptr)->me = NULL;
   return;
 } // end CreateSkyModel

void DeleteSkyModel(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_SkyModel, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Delete: could not recover c struct");
	return;
   }
   ((SkyModel*)ptr)->me = SkyModelUnref(((SkyModel*)ptr)->me);
} // end DeleteSkyModel

ObitSkyModel* SkyModel_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_SkyModel, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Get_me: could not recover c struct");
	return NULL;
   }
   return ((SkyModel*)ptr)->me;
} // end SkyModel_Get_me

void SkyModel_Set_me(PyObject *self, ObitSkyModel *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_SkyModel, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"SkyModel_Set_me: could not recover c struct");
	return;
   }
   ((SkyModel*)ptr)->me =  SkyModelRef(me);
}  // end SkyModel_Set_me


extern ObitSkyModel* SkyModelCopy  (ObitSkyModel *in, ObitSkyModel *out, 
				    ObitErr *err) {
  return ObitSkyModelCopy (in, out, err);
} // end  SkyModelCopy

extern ObitInfoList* SkyModelGetList (ObitSkyModel* in) {
  return ObitInfoListRef(in->info);
}

extern ObitImageMosaic* SkyModelGetImageMosaic (ObitSkyModel* in) {
  return ObitImageMosaicRef(in->mosaic);
}

extern void SkyModelSetImageMosaic (ObitSkyModel* in, ObitImageMosaic *mosaic, 
                              ObitErr *err) {
  in->mosaic = ObitImageMosaicUnref(in->mosaic);  /* Out with the old */
  in->mosaic = ObitImageMosaicRef(mosaic);        /* In with the new */
}

extern ObitSkyModel* SkyModelCreate (char *name, ObitImageMosaic* mosaic) {
 return ObitSkyModelCreate(name, mosaic);
}

extern void SkyModelSubUV (ObitSkyModel* in, ObitUV *inData, ObitUV *outData, ObitErr *err) {
 ObitSkyModelSubUV(in, inData, outData, err);
}

extern void SkyModelDivUV (ObitSkyModel* in, ObitUV *inData, ObitUV *outData, ObitErr *err) {
 ObitSkyModelDivUV(in, inData, outData, err);
}

extern void SkyModelCompressCC (ObitSkyModel* in, ObitErr *err) {
   ObitSkyModelCompressCC (in, err);
}

extern char* SkyModelGetName (ObitSkyModel* in) {
  if (ObitSkyModelIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern long SkyModelIsA (ObitSkyModel* in) {
  return (long)ObitSkyModelIsA(in);
}

%}

