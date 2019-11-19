/* $Id$ */  
/*--------------------------------------------------------------------*/
/* Swig module description for UV data self calibration utilities     */
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
#include "ObitUVSelfCal.h"
#include "ObitUVSoln.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitUVSelfCal *me;
} UVSelfCal;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitUVSelfCal *me;
} UVSelfCal;

%inline %{
extern ObitUVSelfCal* newUVSelfCal (char* name) {
  return newObitUVSelfCal (name);
} // end  newUVSelfCal

extern ObitUVSelfCal* UVSelfCalCreate (char *name, ObitSkyModel *skyModel) {
 return ObitUVSelfCalCreate(name, skyModel);
}

extern ObitUVSelfCal* UVSelfCalUnref (ObitUVSelfCal* in) {
  if (!ObitUVSelfCalIsA(in)) return NULL;
  return ObitUVSelfCalUnref(in);
}

extern ObitUVSelfCal*  UVSelfCalRef (ObitUVSelfCal* in) {
  return ObitUVSelfCalRef(in);
}

void CreateUVSelfCal(PyObject *self, char* name) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UVSelfCal, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Create: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((UVSelfCal*)ptr)->me = newUVSelfCal(name);
   else ((UVSelfCal*)ptr)->me = NULL;
   return;
 } // end CreateUVSelfCal

void DeleteUVSelfCal(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UVSelfCal, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Delete: could not recover c struct");
	return;
   }
   ((UVSelfCal*)ptr)->me = UVSelfCalUnref(((UVSelfCal*)ptr)->me);
} // end DeleteUVSelfCal

ObitUVSelfCal* UVSelfCal_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UVSelfCal, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Get_me: could not recover c struct");
	return NULL;
   }
   return ((UVSelfCal*)ptr)->me;
} // end UVSelfCal_Get_me

void UVSelfCal_Set_me(PyObject *self, ObitUVSelfCal *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UVSelfCal, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"UVSelfCal_Set_me: could not recover c struct");
	return;
   }
   ((UVSelfCal*)ptr)->me =  UVSelfCalRef(me);
}  // end UVSelfCal_Set_me


extern ObitUVSelfCal* UVSelfCalCopy  (ObitUVSelfCal *in, ObitUVSelfCal *out, 
				    ObitErr *err) {
  return ObitUVSelfCalCopy (in, out, err);
} // end  UVSelfCalCopy

extern ObitInfoList* UVSelfCalGetList (ObitUVSelfCal* in) {
  return ObitInfoListRef(in->info);
}

extern ObitSkyModel* UVSelfCalGetSkyModel (ObitUVSelfCal* in) {
  return ObitSkyModelRef(in->skyModel);
}

extern void UVSelfCalSetSkyModel (ObitUVSelfCal* in, ObitSkyModel *skyModel, 
                                     ObitErr *err) {
  in->skyModel = ObitSkyModelUnref(in->skyModel);  /* Out with the old */
  in->skyModel = ObitSkyModelRef(skyModel);        /* In with the new */
}

extern long UVSelfCalIsA (ObitUVSelfCal* in) {
  return (long)ObitUVSelfCalIsA(in);
}

%}


