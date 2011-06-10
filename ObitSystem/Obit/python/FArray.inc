/* $Id$           */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitFarray type                        */
/*                                                                    */
/*;  Copyright (C) 2004-2011                                          */
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

#include "ObitFArray.h"
%}


%inline %{
extern ObitFArray* FArrayCreate(char* name, long ndim, long *naxis) {
   olong i, lnaxis[10];
   for (i=0; i<ndim; i++) lnaxis[i] = (olong)naxis[i];
   return  ObitFArrayCreate (name,(olong) ndim, (olong*)lnaxis);
}

extern float FArrayGetVal(ObitFArray* in, long *pos) {
   float *off, out[4];
   long *Iout;
   olong i, lpos[10];
   for (i=0; i<in->ndim; i++) lpos[i] = (olong)pos[i];
  // check if in bounds
   off = ObitFArrayIndex (in, lpos);
   if (off==NULL) {
	PyErr_SetString(PyExc_RuntimeError,"Position not in array");
        Iout = (long*)out;
        *Iout = ~0;
        return out[0];
   } else {
     out[0] = *off;
   }
   // return NaN rather than Obit magic value
   if (out[0]==ObitMagicF()) {  // create a word with all bits on
        Iout = (long*)out;
        *Iout = ~0;
   }
   return out[0];
}

extern float FArrayGetBlank(void) {
  return ObitMagicF();
}

extern void FArraySetVal(ObitFArray* in, long *pos, float val) {
   float *off;
   olong i, lpos[10];

   for (i=0; i<in->ndim; i++) lpos[i] = (olong)pos[i];

   off = ObitFArrayIndex (in, lpos);
   // check if in bounds
   if (off==NULL) {
      PyErr_SetString(PyExc_RuntimeError,"Position not in array");
      return;
   }
   *off = val;
}

 extern PyObject *FArrayGetBuf (ObitFArray *in) {
   return PyBuffer_FromReadWriteMemory(in->array,
 				      in->arraySize * sizeof(ofloat));
 }
 
extern ObitFArray* FArrayCopy  (ObitFArray *in, ObitFArray *out, ObitErr *err) {
  return ObitFArrayCopy (in, out, err);
} // end FArrayCopy 

extern ObitFArray* FArrayClone  (ObitFArray *in, ObitFArray *out, ObitErr *err) {
  ObitFArrayClone (in, out, err);
  return out;
} // end FArrayCopy 

extern ObitFArray* FArraySubArr  (ObitFArray *in, long *blc, long *trc, 
				  ObitErr *err) {
  olong i, lblc[10], ltrc[10];

  for (i=0; i<in->ndim; i++) lblc[i] = (olong)blc[i];
  for (i=0; i<in->ndim; i++) ltrc[i] = (olong)trc[i];

  return ObitFArraySubArr (in, lblc, ltrc, err);
} // end FArraySubArr

extern ObitFArray* FArrayTranspose  (ObitFArray *in, long *order, ObitErr *err) {
  olong i, lorder[10];
  for (i=0; i<in->ndim; i++) lorder[i] = (olong)order[i];
  return ObitFArrayTranspose (in, lorder, err);
} // end FArrayTranspose

extern int FArrayIsCompatable  (ObitFArray *in1, ObitFArray *in2) {
  return ObitFArrayIsCompatable(in1, in2);
}

extern ObitFArray* FArrayRealloc (ObitFArray* in, long ndim, long *naxis) {
   olong i, lnaxis[10];
   for (i=0; i<ndim; i++) lnaxis[i] = (olong)naxis[i];
   return  ObitFArrayRealloc(in, (olong)ndim, lnaxis);
}

//extern float* FArrayIndex (ObitFArray* in, long *pos) {
//   return ObitFArrayIndex (in, pos);
//}

extern float FArrayMax (ObitFArray* in, long *outValue2) {
   olong i,loutValue2[10];
   float val;
   val = (float)ObitFArrayMax(in, loutValue2);
   for (i=0; i<2; i++) outValue2[i] = (long)loutValue2[i];
   return val;
}

extern float FArrayMaxAbs (ObitFArray* in, long *outValue2) {
   olong i,loutValue2[10];
   float val;
   val = (float)ObitFArrayMaxAbs(in, loutValue2);
   for (i=0; i<2; i++) outValue2[i] = (long)loutValue2[i];
   return val;
}

extern float FArrayMin (ObitFArray* in, long *outValue2) {
   olong i,loutValue2[10];
   float val;
   val = (float)ObitFArrayMin(in, loutValue2);
   for (i=0; i<2; i++) outValue2[i] = (long)loutValue2[i];
   return val;
}

extern void FArrayDeblank (ObitFArray* in, float scalar) {
   ObitFArrayDeblank (in, (ofloat)scalar);
}

extern float FArrayRMS (ObitFArray* in) {
   float out[4];
   long *Iout;
   out[0] = ObitFArrayRMS(in);
   // return NaN rather than Obit magic value
    if (out[0]==ObitMagicF()) {  // create a word with all bits on
        Iout = (long*)&out;
        *Iout = ~0;
   }
   return out[0];
}

extern float FArrayRMS0 (ObitFArray* in) {
   float out[4];
   long *Iout;
   out[0] = ObitFArrayRMS0(in);
   // return NaN rather than Obit magic value
    if (out[0]==ObitMagicF()) {  // create a word with all bits on
        Iout = (long*)&out;
        *Iout = ~0;
   }
   return out[0];
}

extern float FArrayRawRMS (ObitFArray* in) {
   float out[4];
   long *Iout;
   out[0] = ObitFArrayRawRMS(in);
   // return NaN rather than Obit magic value
    if (out[0]==ObitMagicF()) {  // create a word with all bits on
        Iout = (long*)&out;
        *Iout = ~0;
   }
   return out[0];
}

extern float FArrayMode (ObitFArray* in) {
   float out[4];
   long *Iout;
   out[0] = ObitFArrayMode(in);
   // return NaN rather than Obit magic value
   if (out[0]==ObitMagicF()) {  // create a word with all bits on
        Iout = (long*)&out;
        *Iout = ~0;
   }
   return out[0];
}

extern float FArrayMean (ObitFArray* in) {
   float out[4];
   long *Iout;
   out[0] = ObitFArrayMean(in);
   // return NaN rather than Obit magic value
   if (out[0]==ObitMagicF()) {  // create a word with all bits on
        Iout = (long*)&out;
        *Iout = ~0;
   }
   return out[0];
}

extern void FArrayFill (ObitFArray* in, float scalar) {
   return ObitFArrayFill(in, (ofloat)scalar);
}

extern void FArrayNeg (ObitFArray* in) {
   ObitFArrayNeg(in);
}

extern void FArraySin (ObitFArray* in) {
   ObitFArraySin(in);
}

extern void FArrayCos (ObitFArray* in) {
   ObitFArrayCos(in);
}

extern void FArraySqrt (ObitFArray* in) {
   ObitFArraySqrt(in);
}

extern float FArraySum (ObitFArray* in) {
   return ObitFArraySum(in);
}

extern long FArrayCount (ObitFArray* in) {
   return ObitFArrayCount(in);
}

extern void FArraySAdd (ObitFArray* in, float scalar) {
   ObitFArraySAdd(in, (ofloat)scalar);
}

extern void FArraySMul (ObitFArray* in, float scalar) {
   ObitFArraySMul(in, (ofloat)scalar);
}

extern void FArraySDiv (ObitFArray* in, float scalar) {
   ObitFArraySDiv(in, (ofloat)scalar);
}

extern void FArrayClip (ObitFArray* in, float minVal,float maxVal, float newVal) {
   ObitFArrayClip(in, (ofloat)minVal, (ofloat)maxVal, (ofloat)newVal);
}

extern void FArrayInClip (ObitFArray* in, float minVal,float maxVal, float newVal) {
   ObitFArrayInClip(in, (ofloat)minVal, (ofloat)maxVal, (ofloat)newVal);
}

extern void FArrayDivClip (ObitFArray* in1, ObitFArray* in2, float minVal, 
	ObitFArray* out) {
   ObitFArrayDivClip(in1, in2, (ofloat)minVal, out);
}

extern void FArrayClipBlank (ObitFArray* in, float minVal,float maxVal) {
   ofloat fblank = ObitMagicF();
   ObitFArrayClip(in, (ofloat)minVal, (ofloat)maxVal, (ofloat)fblank);
}

extern void FArrayBlank (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArrayBlank (in1, in2, out);
}

extern void FArraySumArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArraySumArr (in1, in2, out);
}

extern void FArrayAvgArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArrayAvgArr (in1, in2, out);
}

extern void FArrayMaxArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArrayMaxArr (in1, in2, out);
}

extern void FArrayMinArr (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArrayMinArr (in1, in2, out);
}

extern void FArrayAdd (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArrayAdd (in1, in2, out);
}

extern void FArraySub (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArraySub (in1, in2, out);
}

extern void FArrayMul (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArrayMul (in1, in2, out);
}

extern void FArrayDiv (ObitFArray* in1, ObitFArray* in2, ObitFArray* out) {
   ObitFArrayDiv (in1, in2, out);
}

extern float FArrayDot (ObitFArray* in1, ObitFArray* in2) {
	return ObitFArrayDot(in1, in2);
}

extern void FArrayMulColRow (ObitFArray* in, ObitFArray* row, ObitFArray* col,
				  ObitFArray* out) {
   ObitFArrayMulColRow (in, row, col, out);
}

extern void FArray2DCenter (ObitFArray* in) {
   ObitFArray2DCenter (in);
}

extern int FArraySymInv2D (ObitFArray* in) {
   int ierr;
   ObitFArray2DSymInv (in, &ierr);
   return ierr;
}

extern void FArrayCGauss2D (ObitFArray* in, long Cen[2], float FWHM) {
   olong lCen[2] = {Cen[0],Cen[1]};
   ObitFArray2DCGauss (in, lCen, (ofloat)FWHM);
}

extern void FArrayEGauss2D (ObitFArray* in, float amp, float Cen[2], float GauMod[3]) {
   ofloat lCen[2] = {Cen[0],Cen[1]};
   ofloat lGauMod[3] = {GauMod[0], GauMod[1], GauMod[2]};
   ObitFArray2DEGauss (in, (float)amp, lCen, lGauMod);
}

extern void FArrayShiftAdd (ObitFArray* in1, long *pos1,
				 ObitFArray* in2, long *pos2, 
				 float scalar, ObitFArray* out) {
   olong i, lpos1[10], lpos2[10];
   for (i=0; i<in1->ndim; i++) lpos1[i] = (olong)pos1[i];
   for (i=0; i<in2->ndim; i++) lpos2[i] = (olong)pos2[i];
   ObitFArrayShiftAdd (in1, lpos1, in2, lpos2, (float)scalar, out);
} // end FArrayShiftAdd

extern void FArrayPad (ObitFArray* in, ObitFArray* out, float factor) {
   ObitFArrayPad (in, out, (float)factor);
} // end FArrayPad

extern void FArraySelInc  (ObitFArray *in, ObitFArray *out, 
	                   long *blc, long *trc, long *inc, ObitErr *err) {
  olong i, lblc[10], ltrc[10], linc[10];

  for (i=0; i<in->ndim; i++) lblc[i] = (olong)blc[i];
  for (i=0; i<in->ndim; i++) ltrc[i] = (olong)trc[i];
  for (i=0; i<in->ndim; i++) linc[i] = (olong)inc[i];
  ObitFArraySelInc (in, out, lblc, ltrc, linc, err);
} // end FArraySelInc

extern ObitFArray* FArrayHisto  (ObitFArray *in, long n, float min, float max) {
  return ObitFArrayHisto (in, (olong)n, (ofloat)min, (ofloat)max);
} // end FArrayHisto

extern ObitInfoList* FArrayGetList (ObitFArray* in) {
  return ObitInfoListRef(in->info);
}

extern char* FArrayGetName (ObitFArray* in) {
  if (!in) return "Undefined";
  return in->name;
} // end  FArrayGetName

extern long FArrayGetNdim (ObitFArray* in) {
  return in->ndim;
} // end  FArrayGetNdim

// returns an array 
extern  PyObject* FArrayGetNaxis (ObitFArray* in) {
  long i;
  PyObject *outList= PyList_New(in->ndim);

  for (i=0; i<in->ndim; i++) {
    PyList_SetItem(outList, i, PyInt_FromLong((long)in->naxis[i]));
  }
  return outList;
} // end  FArrayGetNaxis

extern int FArrayIsA (ObitFArray* in) {
  return ObitFArrayIsA(in);
} // end  FArrayIsA 

ObitFArray* FArrayRef (ObitFArray* in) {
  return ObitFArrayRef (in);
} // end FArrayRef

ObitFArray* FArrayUnref (ObitFArray* in) {
  if (!ObitFArrayIsA(in)) return NULL;
  if (in && (in->ReferenceCount>0)) in = ObitFArrayUnref (in);
  return in;
} // end FArrayUnref

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitFArray *me;
} FArray;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitFArray *me;
} FArray;

%addmethods FArray { 
  FArray(char* name, long ndim, long *naxis) {
     FArray *out;
     out = (FArray *) malloc(sizeof(FArray));
     if (strcmp(name, "None")) out->me = FArrayCreate(name, ndim, naxis);
     else  out->me = NULL;
     return out;
   }
  ~FArray() {
    self->me = FArrayUnref(self->me);
    free(self);
  }
};

