/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for SpectrumFit type                       */
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
#include "ObitSpectrumFit.h"
#include "ObitImage.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitSpectrumFit *me;
} SpectrumFit;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitSpectrumFit *me;
} SpectrumFit;

%inline %{
extern ObitSpectrumFit* newSpectrumFit (char* name) {
  return newObitSpectrumFit (name);
} // end  newSpectrumFit

extern ObitSpectrumFit* SpectrumFitUnref (ObitSpectrumFit* in) {
  if (!ObitSpectrumFitIsA(in)) return NULL;
  return ObitSpectrumFitUnref(in);
}

extern ObitSpectrumFit*  SpectrumFitRef (ObitSpectrumFit* in) {
  return ObitSpectrumFitRef(in);
}

extern ObitSpectrumFit* SpectrumFitCreate (char *name, long nterm) {
  return ObitSpectrumFitCreate((gchar*)name, (olong)nterm);
}

void CreateSpectrumFit(PyObject *self, char* name, long nterm) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_SpectrumFit, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateSpectrumFit: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((SpectrumFit*)ptr)->me = SpectrumFitCreate(name, nterm);
   else  ((SpectrumFit*)ptr)->me = NULL;
   return;
 } // end CreateSpectrumFit

void DeleteSpectrumFit(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_SpectrumFit, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteSpectrumFit: could not recover c struct");
	return;
   }
   ((SpectrumFit*)ptr)->me = ObitSpectrumFitUnref(((SpectrumFit*)ptr)->me);
}// end DeleteSpectrumFit

ObitSpectrumFit* SpectrumFit_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_SpectrumFit, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"SpectrumFit_Get_me: could not recover c struct");
	return NULL;
   }
   return ((SpectrumFit*)ptr)->me;
}  // end SpectrumFit_Get_me

void SpectrumFit_Set_me(PyObject *self, ObitSpectrumFit *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_SpectrumFit, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"SpectrumFit_Set_me: could not recover c struct");
	return;
   }
   ((SpectrumFit*)ptr)->me =  ObitSpectrumFitRef(me);
}  // end SpectrumFit_Set_me

extern ObitSpectrumFit* SpectrumFitCopy  (ObitSpectrumFit *in, ObitSpectrumFit *out, 
				    ObitErr *err) {
  return ObitSpectrumFitCopy (in, out, err);
} // end  SpectrumFitCopy

extern void SpectrumFitCube (ObitSpectrumFit* in, ObitImage *inImage, 
			     ObitImage *outImage, ObitErr *err) {
  ObitSpectrumFitCube(in, inImage, outImage, err);
}

// SWIG REALLY SUCKS: arrays of pointers no longer work
extern void SpectrumFitImArr (ObitSpectrumFit* in, long nimage, 
       ObitImage *im1, ObitImage *im2, ObitImage *im3, ObitImage *im4, ObitImage *im5,
       ObitImage *im6, ObitImage *im7, ObitImage *im8, ObitImage *im9, ObitImage *im10,
       ObitImage *im11, ObitImage *im12, ObitImage *im13, ObitImage *im14, ObitImage *im15,
       ObitImage *im16, ObitImage *im17, ObitImage *im18, ObitImage *im19, ObitImage *im20,
			      ObitImage *outImage, ObitErr *err) {
  ObitImage *imArr[20];
  imArr[0] = im1;
  imArr[1] = im2;
  imArr[2] = im3;
  imArr[3] = im4;
  imArr[4] = im5;
  imArr[5] = im6;
  imArr[6] = im7;
  imArr[7] = im8;
  imArr[8] = im9;
  imArr[9] = im10;
  imArr[10] = im11;
  imArr[11] = im12;
  imArr[12] = im13;
  imArr[13] = im14;
  imArr[14] = im15;
  imArr[15] = im16;
  imArr[16] = im17;
  imArr[17] = im18;
  imArr[18] = im19;
  imArr[19] = im20;
  ObitSpectrumFitImArr(in, (olong)nimage, imArr, outImage, err);
} // end SpectrumFitImArr 

extern void SpectrumFitEval (ObitSpectrumFit* in, ObitImage *inImage,  
                             double outFreq, ObitImage *outImage, 
                             ObitErr *err) {
  ObitSpectrumFitEval(in, inImage, (odouble)outFreq, outImage, err);
}

extern PyObject* SpectrumFitSingle (long nfreq, long nterm, double refFreq, double *freq, float *flux, 
                                    float *sigma, long doBrokePow, ObitErr *err) {
  ofloat *out=NULL;
  olong i, n;
  gboolean ldoBrokePow;
  PyObject *outList=NULL, *o=NULL;

  ldoBrokePow = doBrokePow!=0;
  out = ObitSpectrumFitSingle((olong)nfreq, (olong)nterm, (odouble)refFreq, (odouble*)freq, 
                              (ofloat*)flux, (ofloat*)sigma, ldoBrokePow, err);
  if (err->error) {
        ObitErrLog(err);
        PyErr_SetString(PyExc_TypeError,"Spectral Fit failed");
	o = PyString_FromString("FAILED");
        return o;
  }
  n = 1 + nterm*2;
  outList = PyList_New(n); 
  for (i=0; i<n; i++) {
    o = PyFloat_FromDouble((double)out[i]);
    PyList_SetItem(outList, i, o);
  }
  if (out) g_free(out);
  return outList;
}  // end SpectrumFitPixel

extern ObitInfoList* SpectrumFitGetList (ObitSpectrumFit* in) {
  return ObitInfoListRef(in->info);
}

extern char* SpectrumFitGetName (ObitSpectrumFit* in) {
  if (ObitSpectrumFitIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern long SpectrumFitIsA (ObitSpectrumFit* in) {
  return (long)ObitSpectrumFitIsA(in);
}

%}

