/* $Id:   */  
/*--------------------------------------------------------------------*/
/* Swig module description for SkyModelVMBeam type                    */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

%{
#include "ObitSkyModelVMBeam.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitSkyModelVMBeam *me;
} SkyModelVMBeam;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitSkyModelVMBeam *me;
} SkyModelVMBeam;

%inline %{
extern ObitSkyModelVMBeam* newSkyModelVMBeam (char* name) {
  return newObitSkyModelVMBeam (name);
} // end  newSkyModelVMBeam

extern ObitSkyModelVMBeam* SkyModelVMBeamUnref (ObitSkyModelVMBeam* in) {
  if (!ObitSkyModelVMBeamIsA(in)) return NULL;
  return ObitSkyModelVMBeamUnref(in);
}

extern ObitSkyModelVMBeam*  SkyModelVMBeamRef (ObitSkyModelVMBeam* in) {
  return ObitSkyModelVMBeamRef(in);
}

void CreateSkyModelVMBeam(PyObject *self, char* name) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_SkyModelVMBeam, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Create: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((SkyModelVMBeam*)ptr)->me = newSkyModelVMBeam(name);
   else ((SkyModelVMBeam*)ptr)->me = NULL;
   return;
 } // end CreateSkyModelVMBeam

void DeleteSkyModelVMBeam(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_SkyModelVMBeam, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Delete: could not recover c struct");
	return;
   }
   ((SkyModelVMBeam*)ptr)->me = SkyModelVMBeamUnref(((SkyModelVMBeam*)ptr)->me);
} // end DeleteSkyModelVMBeam

ObitSkyModelVMBeam* SkyModelVMBeam_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_SkyModelVMBeam, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Get_me: could not recover c struct");
	return NULL;
   }
   return ((SkyModelVMBeam*)ptr)->me;
} // end SkyModelVMBeam_Get_me

void SkyModelVMBeam_Set_me(PyObject *self, ObitSkyModelVMBeam *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_SkyModelVMBeam, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"SkyModelVMBeam_Set_me: could not recover c struct");
	return;
   }
   ((SkyModelVMBeam*)ptr)->me =  SkyModelVMBeamRef(me);
}  // end SkyModelVMBeam_Set_me

extern ObitSkyModel* SkyModelVMBeamCast (ObitSkyModelVMBeam *sm) {
  return (ObitSkyModel*)sm;
} // end SkyModelVMBeamCast


extern ObitSkyModelVMBeam* SkyModelVMBeamCopy  (ObitSkyModelVMBeam *in, ObitSkyModelVMBeam *out, 
				    ObitErr *err) {
  return ObitSkyModelVMBeamCopy (in, out, err);
} // end  SkyModelVMBeamCopy

extern ObitInfoList* SkyModelVMBeamGetList (ObitSkyModelVMBeam* in) {
  return ObitInfoListRef(in->info);
}

extern ObitSkyModelVMBeam* SkyModelVMBeamCreate (char *name, ObitImageMosaic* mosaic,
  ObitUV *uvData, ObitImage *IBeam,  ObitImage *VBeam, ObitImage *QBeam,  ObitImage *UBeam, 
  ObitImage *IBeamPh,  ObitImage *VBeamPh, ObitImage *QBeamPh,  ObitImage *UBeamPh, 
  ObitErr *err) {
 return ObitSkyModelVMBeamCreate(name, mosaic, uvData, IBeam, VBeam, QBeam, UBeam, 
    IBeamPh, VBeamPh, QBeamPh, UBeamPh, err);
}

extern char* SkyModelVMBeamGetName (ObitSkyModelVMBeam* in) {
  if (ObitSkyModelVMBeamIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern long SkyModelVMBeamIsA (ObitSkyModelVMBeam* in) {
  return (long)ObitSkyModelVMBeamIsA(in);
}

%}


