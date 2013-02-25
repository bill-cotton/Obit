/* $Id:  $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for RMFit type                             */
/*                                                                    */
/*;  Copyright (C) 2013                                               */
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


%inline %{
extern ObitRMFit* RMFitCopy  (ObitRMFit *in, ObitRMFit *out, 
				    ObitErr *err) {
  return ObitRMFitCopy (in, out, err);
} // end  RMFitCopy

extern ObitRMFit* RMFitUnref (ObitRMFit* in) {
  if (!ObitRMFitIsA(in)) return NULL;
  return ObitRMFitUnref(in);
}

extern ObitRMFit*  RMFitRef (ObitRMFit* in) {
  return ObitRMFitRef(in);
}

extern ObitRMFit* RMFitCreate (char *name, int nterm) {
  return ObitRMFitCreate((gchar*)name, (olong)nterm);
}

extern void RMFitCube (ObitRMFit* in, ObitImage *inQImage, ObitImage *inUImage, 
		       ObitImage *outImage, ObitErr *err) {    
  ObitRMFitCube(in, inQImage, inUImage, outImage, err);
}

extern void RMFitImArr (ObitRMFit* in, int nimage, 
                        ObitImage **imQArr, ObitImage **imUArr, 
			ObitImage *outImage, ObitErr *err) {
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

%addmethods RMFit { 
  RMFit(char* name, int nterm) {
     RMFit *out;
     out = (RMFit *) malloc(sizeof(RMFit));
     if (strcmp(name, "None")) out->me = RMFitCreate(name, nterm);
     else out->me = NULL;
     return out;
   }
  ~RMFit() {
   if (!self) return;  // Not defined
   if (self && self->me && self->me->ReferenceCount>0) {
      self->me = RMFitUnref(self->me);
      free(self);
   }
  }
};

