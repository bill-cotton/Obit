/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitCArray type                        */
/*                                                                    */
/*;  Copyright (C) 2004-2022                                          */
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
#include "ObitCArray.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitCArray *me;
} CArray;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitCArray *me;
} CArray;


%inline %{
extern ObitCArray* CArrayCreate(char* name, long ndim, long *naxis) {
   olong i, lnaxis[10];
   for (i=0; i<ndim; i++) lnaxis[i] = (olong)naxis[i];
   return  ObitCArrayCreate (name, (olong)ndim, lnaxis);
}

ObitCArray* CArrayRef (ObitCArray* in) {
  return ObitCArrayRef (in);
} // end CArrayRef

ObitCArray* CArrayUnref (ObitCArray* in) {
  if (!ObitCArrayIsA(in)) return NULL;
  return ObitCArrayUnref (in);
} // end CArrayUnref

void CreateCArray(PyObject *self, char* name, long ndim, long *longarr) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_CArray, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateCArray: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((CArray*)ptr)->me = CArrayCreate(name, ndim, longarr);
   else  ((CArray*)ptr)->me = NULL;
   return;
 }
void DeleteCArray(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_CArray, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteCArray: could not recover c struct");
	return;
   }
   ((CArray*)ptr)->me = CArrayUnref(((CArray*)ptr)->me);
}// end DeleteCArray

ObitCArray* CArray_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_CArray, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CArray_Get_me: could not recover c struct");
	return NULL;
   }
   return ((CArray*)ptr)->me;
}  // end CArray_Get_me

void CArray_Set_me(PyObject *self, ObitCArray *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_CArray, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CArray_Set_me: could not recover c struct");
	return;
   }
   ((CArray*)ptr)->me =  ObitCArrayRef(me);
}  // end CArray_Set_me

extern ObitCArray* CArrayCopy  (ObitCArray *in, ObitCArray *out, ObitErr *err) {
  return ObitCArrayCopy (in, out, err);
} // end CArrayCopy 

extern int CArrayIsCompatable  (ObitCArray *in1, ObitCArray *in2) {
  return ObitCArrayIsCompatable(in1, in2);
}

extern ObitCArray* CArrayRealloc (ObitCArray* in, long ndim, long *naxis) {
   olong i, lnaxis[10];
   for (i=0; i<ndim; i++) lnaxis[i] = (olong)naxis[i];
   return  ObitCArrayRealloc(in, (olong)ndim, lnaxis);
}

extern PyObject* CArrayGetVal(ObitCArray* in, long *pos, float *val) {
   float *off, out[2];
   olong i, lpos[10];
   long *Iout;
   PyObject *outList, *o;

   for (i=0; i<in->ndim; i++) lpos[i] = (olong)pos[i];
  // check if in bounds
   off = ObitCArrayIndex (in, lpos);
   if (off==NULL) {
	PyErr_SetString(PyExc_RuntimeError,"Position not in array");
        Iout = (long*)out;  /* I think this gives a NaN */
        *Iout = ~0;
        out[1] = out[0];
   } else {
     out[0] = off[0];out[1] = off[1];
   }
   val[0] = out[0];
   val[1] = out[1];
   outList = PyList_New(2);
   o = PyFloat_FromDouble((double)out[0]);
   PyList_SetItem(outList, 0, o);
   o = PyFloat_FromDouble((double)out[1]);
   PyList_SetItem(outList, 1, o);	
   return outList;
}

extern void CArraySetVal(ObitCArray* in, long *pos, float *val) {
   float *off;
   olong i, lpos[10];

   for (i=0; i<in->ndim; i++) lpos[i] = (olong)pos[i];
   off = ObitCArrayIndex (in, lpos);
   // check if in bounds
   if (off==NULL) {
      PyErr_SetString(PyExc_RuntimeError,"Position not in array");
      return;
   }
   off[0] = val[0];
   off[1] = val[1];
}

extern float CArrayMaxAbs (ObitCArray* in, long *outValue2) {
   olong i,loutValue2[10];
   float val;
   val = (float)ObitCArrayMaxAbs(in, loutValue2);
   for (i=0; i<2; i++) outValue2[i] = (long)loutValue2[i];
	
   return val;
}

extern void CArrayNeg (ObitCArray* in) {
   ObitCArrayNeg(in);
}

extern void CArrayConjg (ObitCArray* in) {
   ObitCArrayConjg(in);
}

extern void CArrayFill (ObitCArray* in, float *cmpx) {
   ObitCArrayFill(in, cmpx);
}

extern void CArraySAdd (ObitCArray* in, float scalar) {
   ObitCArraySAdd(in, scalar);
}

extern void CArraySMul (ObitCArray* in, float scalar) {
   ObitCArraySMul(in, scalar);
}

extern void CArrayCSAdd (ObitCArray* in, float *scalar) {
   ObitCArrayCSAdd(in, scalar);
}

extern void CArrayCSMul (ObitCArray* in, float *scalar) {
   ObitCArrayCSMul(in, scalar);
}

extern void CArrayAdd (ObitCArray* in1, ObitCArray* in2, ObitCArray* out) {
   ObitCArrayAdd (in1, in2, out);
}

extern void CArraySub (ObitCArray* in1, ObitCArray* in2, ObitCArray* out) {
   ObitCArraySub (in1, in2, out);
}

extern void CArrayMul (ObitCArray* in1, ObitCArray* in2, ObitCArray* out) {
   ObitCArrayMul (in1, in2, out);
}

extern void CArrayDiv (ObitCArray* in1, ObitCArray* in2, ObitCArray* out) {
   ObitCArrayDiv (in1, in2, out);
}

extern ObitFArray* CArrayMakeF  (ObitCArray *Cin) {
  return ObitCArrayMakeF  (Cin);
}

extern ObitCArray* CArrayMakeC  (ObitFArray *Fin) {
  return ObitCArrayMakeC  (Fin);
}

extern int CArrayIsFCompatable  (ObitCArray *Cin, ObitFArray *Fin) {
  return ObitCArrayIsFCompatable  (Cin, Fin);
}

extern void CArrayFMul (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out) {
  ObitCArrayFMul (Cin, Fin, out);
}

extern void CArrayFDiv (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out) {
  ObitCArrayFDiv (Cin, Fin, out);
}

extern void CArrayFAdd (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out) {
  ObitCArrayFAdd (Cin, Fin, out);
}

extern void CArrayFRot (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out) {
  ObitCArrayFRot (Cin, Fin, out);
}

extern void CArrayComplex (ObitFArray* Fin1, ObitFArray* Fin2, ObitCArray* out) {
  ObitCArrayComplex (Fin1, Fin2, out);
}

extern void CArrayReal (ObitCArray* in, ObitFArray* out) {
  ObitCArrayReal (in, out);
}

extern void CArrayImag (ObitCArray* in, ObitFArray* out) {
  ObitCArrayImag (in, out);
}

extern void CArrayAmp (ObitCArray* in, ObitFArray* out) {
  ObitCArrayAmp (in, out);
}

extern void CArrayPhase (ObitCArray* in, ObitFArray* out) {
  ObitCArrayPhase (in, out);
}

extern void CArray2DCenter(ObitCArray* in) {
  ObitCArray2DCenter (in);
}

extern ObitCArray* CArrayAddConjg (ObitCArray* in, long numConjRow) {
  return ObitCArrayAddConjg (in, (olong)numConjRow);
}

extern char* CArrayGetName (ObitCArray* in) {
  return in->name;
} // end  CArrayGetName

extern long CArrayGetNdim (ObitCArray* in) {
  return in->ndim;
} // end  CArrayGetNdim

// returns an array of 7 elements no matter
extern PyObject* CArrayGetNaxis (ObitCArray* in) {
  long i;
  PyObject *outList= PyList_New(in->ndim);

  for (i=0; i<in->ndim; i++) {
    PyList_SetItem(outList, i, PyLong_FromLong((long)in->naxis[i]));
  }
  return outList;
  return outList;
} // end  CArrayGetNaxis

extern int CArrayIsA (ObitCArray* in) {
  return ObitCArrayIsA(in);
} // end  CArrayIsA 

%}


