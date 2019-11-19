/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for BeamShape type                         */
/*                                                                    */
/*;  Copyright (C) 2008-2019                                          */
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
#include "ObitBeamShape.h"
#include "ObitImage.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitBeamShape *me;
} BeamShape;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitBeamShape *me;
} BeamShape;


%inline %{
extern ObitBeamShape* newBeamShape (char* name) {
  return newObitBeamShape (name);
} // end  newBeamShape

extern ObitBeamShape* BeamShapeUnref (ObitBeamShape* in) {
  if (!ObitBeamShapeIsA(in)) return NULL;
  return ObitBeamShapeUnref(in);
}

extern ObitBeamShape*  BeamShapeRef (ObitBeamShape* in) {
  return ObitBeamShapeRef(in);
}

extern ObitBeamShape* BeamShapeCreate (char *name, ObitImage *image, 
				       float pbmin, float antSize, 
				       int doGain) {
  gboolean ldoGain = (doGain!=0);
  return ObitBeamShapeCreate((gchar*)name, image, (ofloat)pbmin, (ofloat)antSize,
			     ldoGain);
} // end BeamShapeCreate

void CreateBeamShape(PyObject *self, char* name, ObitImage *image, float pbmin, float antSize,
            long doGain) {
   void *ptr;
   gboolean ldoGain = (doGain!=0);
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_BeamShape, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateBeamShape: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((BeamShape*)ptr)->me = BeamShapeCreate(name, image,
                                                         (ofloat)pbmin, (ofloat)antSize,
                                			 ldoGain);
   else  ((BeamShape*)ptr)->me = NULL;
   return;
 } // end CreateBeamShape

void DeleteBeamShape(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_BeamShape, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteBeamShape: could not recover c struct");
	return;
   }
   ((BeamShape*)ptr)->me = ObitBeamShapeUnref(((BeamShape*)ptr)->me);
}// end DeleteBeamShape

ObitBeamShape* BeamShape_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_BeamShape, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"BeamShape_Get_me: could not recover c struct");
	return NULL;
   }
   return ((BeamShape*)ptr)->me;
}  // end BeamShape_Get_me

void BeamShape_Set_me(PyObject *self, ObitBeamShape *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_BeamShape, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"BeamShape_Set_me: could not recover c struct");
	return;
   }
   ((BeamShape*)ptr)->me =  ObitBeamShapeRef(me);
}  // end BeamShape_Set_me

extern ObitBeamShape* BeamShapeCopy  (ObitBeamShape *in, ObitBeamShape *out, 
				    ObitErr *err) {
  return ObitBeamShapeCopy (in, out, err);
} // end  BeamShapeCopy

extern float BeamShapeGain (ObitBeamShape* in, double ra, double dec,
		            float parAng ) {
  return ObitBeamShapeGain(in, (odouble)ra, (odouble)dec, (ofloat)parAng);
} // end BeamShapeGain

extern float BeamShapeGainSym (ObitBeamShape* in, double Angle) {
  return ObitBeamShapeGainSym(in, (odouble)Angle);
} // end BeamShapeGainSym

extern char* BeamShapeGetName (ObitBeamShape* in) {
  if (ObitBeamShapeIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int BeamShapeIsA (ObitBeamShape* in) {
  return ObitBeamShapeIsA(in);
}

extern float BeamShapeAngle (ObitBeamShape* in, float ra, float dec, float parAng) {
  return (float) ObitBeamShapeAngle (in, (odouble)ra, (odouble)dec, (ofloat)parAng);
}

extern float BeamShapeGetFreq (ObitBeamShape* in) {
  return (float) ObitBeamShapeGetFreq (in);
}

extern void BeamShapeSetFreq (ObitBeamShape* in, float newFreq) {
  ObitBeamShapeSetFreq (in, (odouble)newFreq);
}

%}

