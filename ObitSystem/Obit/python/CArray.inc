/* $Id: CArray.inc,v 1.7 2006/10/11 17:24:48 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitCArray type                        */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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


%inline %{
extern ObitCArray* CArrayCreate(char* name, long ndim, long *naxis) {
   olong i, lnaxis[10];
   for (i=0; i<ndim; i++) lnaxis[i] = (olong)naxis[i];
   return  ObitCArrayCreate (name, (olong)ndim, lnaxis);
}

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

extern void CArrayGetVal(ObitCArray* in, long *pos, float val[2]) {
   float *out;
   olong i, lpos[10];

   for (i=0; i<in->ndim; i++) lpos[i] = (olong)pos[i];
   out = ObitCArrayIndex(in, lpos);
   val[0] = out[0];
   val[1] = out[1];
}

extern void CArraySetVal(ObitCArray* in, long *pos, float val[2]) {
   float *off;
   olong i, lpos[10];

   for (i=0; i<in->ndim; i++) lpos[i] = (olong)pos[i];
   off = ObitCArrayIndex (in, lpos);
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

extern void CArrayFill (ObitCArray* in, float cmpx[2]) {
   ObitCArrayFill(in, cmpx);
}

extern void CArraySAdd (ObitCArray* in, float scalar) {
   ObitCArraySAdd(in, scalar);
}

extern void CArraySMul (ObitCArray* in, float scalar) {
   ObitCArraySMul(in, scalar);
}

extern void CArrayCSAdd (ObitCArray* in, float scalar[2]) {
   ObitCArrayCSAdd(in, scalar);
}

extern void CArrayCSMul (ObitCArray* in, float scalar[2]) {
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
    PyList_SetItem(outList, i, PyInt_FromLong((long)in->naxis[i]));
  }
  return outList;
  return outList;
} // end  CArrayGetNaxis

extern int CArrayIsA (ObitCArray* in) {
  return ObitCArrayIsA(in);
} // end  CArrayIsA 

ObitCArray* CArrayRef (ObitCArray* in) {
  return ObitCArrayRef (in);
} // end CArrayRef

ObitCArray* CArrayUnref (ObitCArray* in) {
  if (!ObitCArrayIsA(in)) return NULL;
  return ObitCArrayUnref (in);
} // end CArrayUnref



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

%addmethods CArray { 
  CArray(char* name, long ndim, long *naxis) {
     CArray *out;
     out = (CArray *) malloc(sizeof(CArray));
     if (strcmp(name, "None")) out->me = CArrayCreate(name, ndim, naxis);
     else out->me = NULL;
     return out;
   }
  ~CArray() {
    self->me = CArrayUnref(self->me);
    free(self);
  }
};

