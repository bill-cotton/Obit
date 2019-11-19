/* $Id$  */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitFFT type                           */
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
#include "ObitFFT.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitFFT *me;
} FFT;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitFFT *me;
} FFT;


%inline %{
// dir   1 => OBIT_FFT_Forward(R2C), 2 => OBIT_FFT_Reverse(C2R)
// type  1 => OBIT_FFT_FullComplex,  2 => OBIT_FFT_HalfComplex
extern ObitFFT* FFTCreate(char* name, long dir, long type, long rank, long *dim) {
   olong i, ldim[7];
   ObitFFTdir ldir;
   ObitFFTtype ltype;
   if (dir==1) ldir = OBIT_FFT_Forward;
   else ldir = OBIT_FFT_Reverse;
   if (type==1) ltype = OBIT_FFT_FullComplex;
   else ltype = OBIT_FFT_HalfComplex;
   for (i=0; i<rank; i++) ldim[i] = (olong)dim[i];
   return newObitFFT  (name, ldir, ltype, (olong)rank, ldim);
}

ObitFFT* FFTRef (ObitFFT* in) {
  return ObitFFTRef (in);
} // end FFTRef

ObitFFT* FFTUnref (ObitFFT* in) {
  if (!ObitFFTIsA(in)) return NULL;
  return ObitFFTUnref (in);
} // end FFTUnref

void CreateFFT(PyObject *self, char* name, long dir, long type, long rank, long *dim) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_FFT, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateFFT: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((FFT*)ptr)->me = FFTCreate(name, dir, type, rank, dim);
   else  ((FFT*)ptr)->me = NULL;
   return;
 } // end CreateFFT

void DeleteFFT(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_FFT, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteFFT: could not recover c struct");
	return;
   }
   ((FFT*)ptr)->me = ObitFFTUnref(((FFT*)ptr)->me);
}// end DeleteFFT

ObitFFT* FFT_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_FFT, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"FFT_Get_me: could not recover c struct");
	return NULL;
   }
   return ((FFT*)ptr)->me;
}  // end FFT_Get_me

void FFT_Set_me(PyObject *self, ObitFFT *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_FFT, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"FFT_Set_me: could not recover c struct");
	return;
   }
   ((FFT*)ptr)->me =  ObitFFTRef(me);
}  // end FFT_Set_me

extern long FFTSuggestSize  (long length) {
  return (long)ObitFFTSuggestSize(length);
}

extern void FFTR2C(ObitFFT* in, ObitFArray *inArray, ObitCArray *outArray) {
   ObitFFTR2C(in, inArray, outArray);
}

extern void FFTC2R(ObitFFT* in, ObitCArray *inArray, ObitFArray *outArray) {
   ObitFFTC2R(in, inArray, outArray);
}

extern void FFTC2C(ObitFFT* in, ObitCArray *inArray, ObitCArray *outArray) {
   ObitFFTC2C(in, inArray, outArray);
}

extern char* FFTGetName (ObitFFT* in) {
  return in->name;
} // end  FFTGetName

extern long FFTGetRank (ObitFFT* in) {
  return (long)in->rank;
} // end  FFTGetRank

// returns a list of longs
extern PyObject* FFTGetDim (ObitFFT* in) {
  olong i;
  PyObject *olist = PyList_New(IM_MAXDIM);
  for (i=0; i<in->rank; i++)
     PyList_SetItem(olist, i, PyLong_FromLong((long)in->dim[i]));
  return olist;
} // end  FFTGetDim

extern long FFTIsA (ObitFFT* in) {
  return (olong)ObitFFTIsA(in);
} // end  FFTIsA 



%}

