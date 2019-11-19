/* $Id: $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for FullBeam type                          */
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
#include "ObitFullBeam.h"
#include "ObitImage.h"
#include "ObitImageDesc.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitFullBeam *me;
} FullBeam;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitFullBeam *me;
} FullBeam;

%inline %{
extern ObitFullBeam* newFullBeam (char* name) {
  return newObitFullBeam (name);
} // end  newFullBeam

extern ObitFullBeam* FullBeamCreate (char *name, ObitImage *image, 
				     ObitErr *err) {
  return ObitFullBeamCreate((gchar*)name, NULL, image, err);
} // end FullBeamCreate

extern ObitFullBeam* FullBeamUnref (ObitFullBeam* in) {
  if (!ObitFullBeamIsA(in)) return NULL;
  return ObitFullBeamUnref(in);
}

extern ObitFullBeam*  FullBeamRef (ObitFullBeam* in) {
  return ObitFullBeamRef(in);
}

void CreateFullBeam(PyObject *self, char* name, ObitImage *image, ObitErr *err) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_FullBeam, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateFullBeam: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((FullBeam*)ptr)->me = FullBeamCreate((gchar*)name, image, err);
   else  ((FullBeam*)ptr)->me = NULL;
   return;
 } // end CreateFullBeam

void DeleteFullBeam(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_FullBeam, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteFullBeam: could not recover c struct");
	return;
   }
   ((FullBeam*)ptr)->me = ObitFullBeamUnref(((FullBeam*)ptr)->me);
}// end DeleteFullBeam

ObitFullBeam* FullBeam_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_FullBeam, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"FullBeam_Get_me: could not recover c struct");
	return NULL;
   }
   return ((FullBeam*)ptr)->me;
}  // end FullBeam_Get_me

void FullBeam_Set_me(PyObject *self, ObitFullBeam *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_FullBeam, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"FullBeam_Set_me: could not recover c struct");
	return;
   }
   ((FullBeam*)ptr)->me =  ObitFullBeamRef(me);
}  // end FullBeam_Set_me

extern ObitFullBeam* FullBeamCopy  (ObitFullBeam *in, ObitFullBeam *out, 
				    ObitErr *err) {
  return ObitFullBeamCopy (in, out, err);
} // end  FullBeamCopy

extern float FullBeamValue (ObitFullBeam* in, double dra, double ddec,
		            float parAng, long plane, ObitErr *err) {
  return ObitFullBeamValue(in, (odouble)dra, (odouble)ddec, (ofloat)parAng, (olong)plane, err);
} // end FullBeamGain

extern long FullBeamFindPlane (ObitFullBeam* in, double freq) {
  return (long)ObitFullBeamFindPlane(in, (odouble)freq);
} // end FullBeamFindPlane

extern char* FullBeamGetName (ObitFullBeam* in) {
  if (ObitFullBeamIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern ObitInfoList* FullBeamGetList (ObitFullBeam* in) {
  return ObitInfoListRef(in->info);
}

extern ObitImageDesc* FullBeamGetDesc (ObitFullBeam* in) {
  return ObitImageDescRef(in->BeamDesc);
}

extern long FullBeamIsA (ObitFullBeam* in) {
  return (long)ObitFullBeamIsA(in);
}

%}

