/* $Id: SpectrumFit.inc,v 1.4 2008/05/16 23:29:02 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for SpectrumFit type                       */
/*                                                                    */
/*;  Copyright (C) 2008                                               */
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


%inline %{
extern ObitSpectrumFit* newSpectrumFit (char* name) {
  return newObitSpectrumFit (name);
} // end  newSpectrumFit

extern ObitSpectrumFit* SpectrumFitCopy  (ObitSpectrumFit *in, ObitSpectrumFit *out, 
				    ObitErr *err) {
  return ObitSpectrumFitCopy (in, out, err);
} // end  SpectrumFitCopy

extern ObitSpectrumFit* SpectrumFitUnref (ObitSpectrumFit* in) {
  if (!ObitSpectrumFitIsA(in)) return NULL;
  return ObitSpectrumFitUnref(in);
}

extern ObitSpectrumFit*  SpectrumFitRef (ObitSpectrumFit* in) {
  return ObitSpectrumFitRef(in);
}

extern ObitSpectrumFit* SpectrumFitCreate (char *name, int nterm) {
  return ObitSpectrumFitCreate((gchar*)name, (olong)nterm);
}

extern void SpectrumFitCube (ObitSpectrumFit* in, ObitImage *inImage, 
			     ObitImage *outImage, ObitErr *err) {
  ObitSpectrumFitCube(in, inImage, outImage, err);
}

extern void SpectrumFitImArr (ObitSpectrumFit* in, int nimage, ObitImage **imArr, 
			      ObitImage *outImage, ObitErr *err) {
  ObitSpectrumFitImArr(in, (olong)nimage, imArr, outImage, err);
} // end SpectrumFitImArr 

extern void SpectrumFitEval (ObitSpectrumFit* in, ObitImage *inImage,  
                             double outFreq, ObitImage *outImage, 
                             ObitErr *err) {
  ObitSpectrumFitEval(in, inImage, (odouble)outFreq, outImage, err);
}

extern PyObject* SpectrumFitSingle (int nfreq, int nterm, double *freq, float *flux, float *sigma,
                                 ObitErr *err) {
  ofloat *out=NULL;
  olong i, n;
  PyObject *outList=NULL, *o=NULL;

  out = ObitSpectrumFitSingle((olong)nfreq, (olong)nterm, (odouble*)freq, (ofloat*)flux, (ofloat*)sigma, err);
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

extern int SpectrumFitIsA (ObitSpectrumFit* in) {
  return ObitSpectrumFitIsA(in);
}

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

%addmethods SpectrumFit { 
  SpectrumFit(char* name, int nterm) {
     SpectrumFit *out;
     out = (SpectrumFit *) malloc(sizeof(SpectrumFit));
     if (strcmp(name, "None")) out->me = SpectrumFitCreate(name, nterm);
     else out->me = NULL;
     return out;
   }
  ~SpectrumFit() {
   if (!self) return;  // Not defined
   if (self && self->me && self->me->ReferenceCount>0) {
      self->me = SpectrumFitUnref(self->me);
      free(self);
   }
  }
};

