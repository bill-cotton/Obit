/* $Id$ */  
/*--------------------------------------------------------------------*/
/* Swig module description for UV data self calibration utilities     */
/*                                                                    */
/*;  Copyright (C) 2006,2019                                          */
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
#include "ObitUVGSolve.h"
#include "ObitUVSoln.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitUVGSolve *me;
} UVGSolve;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitUVGSolve *me;
} UVGSolve;


%inline %{
extern ObitUVGSolve* newUVGSolve (char* name) {
  return newObitUVGSolve (name);
} // end  newUVGSolve

extern ObitUVGSolve* UVGSolveCreate (char *name) {
 return ObitUVGSolveCreate(name);
}

extern ObitUVGSolve* UVGSolveUnref (ObitUVGSolve* in) {
  if (!ObitUVGSolveIsA(in)) return NULL;
  return ObitUVGSolveUnref(in);
}

extern ObitUVGSolve*  UVGSolveRef (ObitUVGSolve* in) {
  return ObitUVGSolveRef(in);
}

void CreateUVGSolve(PyObject *self, char* name) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UVGSolve, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Create: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((UVGSolve*)ptr)->me = newUVGSolve(name);
   else ((UVGSolve*)ptr)->me = NULL;
   return;
 } // end CreateUVGSolve

void DeleteUVGSolve(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UVGSolve, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Delete: could not recover c struct");
	return;
   }
   ((UVGSolve*)ptr)->me = UVGSolveUnref(((UVGSolve*)ptr)->me);
} // end DeleteUVGSolve

ObitUVGSolve* UVGSolve_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UVGSolve, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Get_me: could not recover c struct");
	return NULL;
   }
   return ((UVGSolve*)ptr)->me;
} // end UVGSolve_Get_me

void UVGSolve_Set_me(PyObject *self, ObitUVGSolve *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_UVGSolve, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"UVGSolve_Set_me: could not recover c struct");
	return;
   }
   ((UVGSolve*)ptr)->me =  UVGSolveRef(me);
}  // end TimeFilter_Set_me

extern ObitUVGSolve* UVGSolveCopy  (ObitUVGSolve *in, ObitUVGSolve *out, 
				    ObitErr *err) {
  return ObitUVGSolveCopy (in, out, err);
} // end  UVGSolveCopy

extern ObitInfoList* UVGSolveGetList (ObitUVGSolve* in) {
  return ObitInfoListRef(in->info);
}

extern long UVGSolveIsA (ObitUVGSolve* in) {
  return (long)ObitUVGSolveIsA(in);
}

extern ObitTable* UVGSolveCal (ObitUVGSolve *in, ObitUV *inUV, ObitUV *outUV, 
                               ObitUVSel *sel, ObitErr *err)
{
  return (ObitTable*)ObitUVGSolveCal (in, inUV, outUV, sel, err);
} // end UVGSolveCal


%}

