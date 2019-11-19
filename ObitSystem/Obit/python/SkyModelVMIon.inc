/* $Id:   */  
/*--------------------------------------------------------------------*/
/* Swig module description for SkyModelVMIon type                     */
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
#include "ObitSkyModelVMIon.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitSkyModelVMIon *me;
} SkyModelVMIon;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitSkyModelVMIon *me;
} SkyModelVMIon;


%inline %{
extern ObitSkyModelVMIon* newSkyModelVMIon (char* name) {
  return newObitSkyModelVMIon (name);
} // end  newSkyModelVMIon

extern ObitSkyModelVMIon* SkyModelVMIonUnref (ObitSkyModelVMIon* in) {
  if (!ObitSkyModelVMIonIsA(in)) return NULL;
  return ObitSkyModelVMIonUnref(in);
}

extern ObitSkyModelVMIon*  SkyModelVMIonRef (ObitSkyModelVMIon* in) {
  return ObitSkyModelVMIonRef(in);
}

void CreateSkyModelVMIon(PyObject *self, char* name) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_SkyModelVMIon, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Create: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((SkyModelVMIon*)ptr)->me = newSkyModelVMIon(name);
   else ((SkyModelVMIon*)ptr)->me = NULL;
   return;
 } // end CreateSkyModelVMIon

void DeleteSkyModelVMIon(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_SkyModelVMIon, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Delete: could not recover c struct");
	return;
   }
   ((SkyModelVMIon*)ptr)->me = SkyModelVMIonUnref(((SkyModelVMIon*)ptr)->me);
} // end DeleteSkyModelVMIon

ObitSkyModelVMIon* SkyModelVMIon_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_SkyModelVMIon, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Get_me: could not recover c struct");
	return NULL;
   }
   return ((SkyModelVMIon*)ptr)->me;
} // end SkyModelVMIon_Get_me

void SkyModelVMIon_Set_me(PyObject *self, ObitSkyModelVMIon *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_SkyModelVMIon, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"SkyModelVMIon_Set_me: could not recover c struct");
	return;
   }
   ((SkyModelVMIon*)ptr)->me =  SkyModelVMIonRef(me);
}  // end SkyModelVMIon_Set_me

extern ObitSkyModel* SkyModelVMIonCast (ObitSkyModelVMIon *sm) {
  return (ObitSkyModel*)sm;
} // end SkyModelVMIonCast

extern ObitSkyModelVMIon* SkyModelVMIonCopy  (ObitSkyModelVMIon *in, ObitSkyModelVMIon *out, 
				    ObitErr *err) {
  return ObitSkyModelVMIonCopy (in, out, err);
} // end  SkyModelVMIonCopy

extern ObitInfoList* SkyModelVMIonGetList (ObitSkyModelVMIon* in) {
  return ObitInfoListRef(in->info);
}

extern ObitSkyModelVMIon* SkyModelVMIonCreate (char *name, ObitImageMosaic* mosaic) {
 return ObitSkyModelVMIonCreate(name, mosaic);
}

extern char* SkyModelVMIonGetName (ObitSkyModelVMIon* in) {
  if (ObitSkyModelVMIonIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern long SkyModelVMIonIsA (ObitSkyModelVMIon* in) {
  return (long)ObitSkyModelVMIonIsA(in);
}

%}


