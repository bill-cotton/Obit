/* $Id:  $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for RMFit type                             */
/*                                                                    */
/*;  Copyright (C) 2013,2019                                          */
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
#include "ObitRMFit.h"
#include "ObitImage.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitRMFit *me;
} RMFit;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitRMFit *me;
} RMFit;


%inline %{
extern ObitRMFit* RMFitUnref (ObitRMFit* in) {
  if (!ObitRMFitIsA(in)) return NULL;
  return ObitRMFitUnref(in);
}

extern ObitRMFit*  RMFitRef (ObitRMFit* in) {
  return ObitRMFitRef(in);
}

extern ObitRMFit* RMFitCreate (char *name, long nterm) {
  return ObitRMFitCreate((gchar*)name, (olong)nterm);
}

extern ObitRMFit* RMFitCopy  (ObitRMFit *in, ObitRMFit *out, 
				    ObitErr *err) {
  return ObitRMFitCopy (in, out, err);
} // end  RMFitCopy

void CreateRMFit(PyObject *self, char* name, long nterm) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_RMFit, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateRMFit: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((RMFit*)ptr)->me = RMFitCreate(name, nterm);
   else  ((RMFit*)ptr)->me = NULL;
   return;
 } // end CreateRMFit

void DeleteRMFit(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_RMFit, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteRMFit: could not recover c struct");
	return;
   }
   ((RMFit*)ptr)->me = ObitRMFitUnref(((RMFit*)ptr)->me);
}// end DeleteRMFit

ObitRMFit* RMFit_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_RMFit, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"RMFit_Get_me: could not recover c struct");
	return NULL;
   }
   return ((RMFit*)ptr)->me;
}  // end RMFit_Get_me

void RMFit_Set_me(PyObject *self, ObitRMFit *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_RMFit, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"RMFit_Set_me: could not recover c struct");
	return;
   }
   ((RMFit*)ptr)->me =  ObitRMFitRef(me);
}  // end RMFit_Set_me


extern void RMFitCube (ObitRMFit* in, ObitImage *inQImage, ObitImage *inUImage, 
		       ObitImage *outImage, ObitErr *err) {    
  ObitRMFitCube(in, inQImage, inUImage, outImage, err);
}

// SWIG REALLY SUCKS: arrays of pointers no longer work
extern void RMFitImArr (ObitRMFit* in, int nimage, 
       ObitImage *imQ1, ObitImage *imQ2, ObitImage *imQ3, ObitImage *imQ4, ObitImage *imQ5,
       ObitImage *imQ6, ObitImage *imQ7, ObitImage *imQ8, ObitImage *imQ9, ObitImage *imQ10,
       ObitImage *imQ11, ObitImage *imQ12, ObitImage *imQ13, ObitImage *imQ14, ObitImage *imQ15,
       ObitImage *imQ16, ObitImage *imQ17, ObitImage *imQ18, ObitImage *imQ19, ObitImage *imQ20,
       ObitImage *imU1, ObitImage *imU2, ObitImage *imU3, ObitImage *imU4, ObitImage *imU5,
       ObitImage *imU6, ObitImage *imU7, ObitImage *imU8, ObitImage *imU9, ObitImage *imU10,
       ObitImage *imU11, ObitImage *imU12, ObitImage *imU13, ObitImage *imU14, ObitImage *imU15,
       ObitImage *imU16, ObitImage *imU17, ObitImage *imU18, ObitImage *imU19, ObitImage *imU20,
			ObitImage *outImage, ObitErr *err) {
  ObitImage *imQArr[20], *imUArr[20];
  imQArr[0] = imQ1; imUArr[0] = imU1;
  imQArr[1] = imQ2; imUArr[1] = imU2;
  imQArr[2] = imQ3; imUArr[2] = imU3;
  imQArr[3] = imQ4; imUArr[3] = imU4;
  imQArr[4] = imQ5; imUArr[4] = imU5;
  imQArr[5] = imQ6; imUArr[5] = imU6;
  imQArr[6] = imQ7; imUArr[6] = imU7;
  imQArr[7] = imQ8; imUArr[7] = imU8;
  imQArr[8] = imQ9; imUArr[8] = imU9;
  imQArr[9] = imQ10; imUArr[9] = imU10;
  imQArr[10] = imQ11; imUArr[10] = imU11;
  imQArr[11] = imQ12; imUArr[11] = imU12;
  imQArr[12] = imQ13; imUArr[12] = imU13;
  imQArr[13] = imQ14; imUArr[13] = imU14;
  imQArr[14] = imQ15; imUArr[14] = imU15;
  imQArr[15] = imQ16; imUArr[15] = imU16;
  imQArr[16] = imQ17; imUArr[16] = imU17;
  imQArr[17] = imQ18; imUArr[17] = imU18;
  imQArr[18] = imQ19; imUArr[18] = imU19;
  imQArr[19] = imQ20; imUArr[19] = imU20;
  ObitRMFitImArr(in, (olong)nimage, imQArr, imUArr, outImage, err);
} // end RMFitImArr 
 

extern ObitRMFit* newRMFit (char* name) {
  return newObitRMFit (name);
} // end  newRMFit

extern ObitInfoList* RMFitGetList (ObitRMFit* in) {
  return ObitInfoListRef(in->info);
}

extern char* RMFitGetName (ObitRMFit* in) {
  if (ObitRMFitIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int RMFitIsA (ObitRMFit* in) {
  return ObitRMFitIsA(in);
}

extern PyObject* RMFitSingle (int nfreq, int nterm, double refLamb2, double *lamb2, 
                              float *qflux, float *qsigma, float *uflux, float *usigma, 
                              ObitErr *err) {
  ofloat *out=NULL;
  olong i, n;
  PyObject *outList=NULL, *o=NULL;
  out = ObitRMFitSingle((olong)nfreq, (olong)nterm, (odouble)refLamb2, (odouble*)lamb2, 
                              (ofloat*)qflux, (ofloat*)qsigma, (ofloat*)uflux, (ofloat*)usigma, err);
  if (err->error) {
        ObitErrLog(err);
        PyErr_SetString(PyExc_TypeError,"RM Fit failed");
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
}  // end RMFitSingle

%}

