/* $Id: $           */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitGPUFarray type                     */
/*                                                                    */
/*;  Copyright (C) 2014                                               */
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

#include "ObitGPUFArray.h"
%}


%inline %{
extern ObitGPUFArray* GPUFArrayCreate(char* name, long ndim, long *naxis) {
   olong i, lnaxis[10];
   for (i=0; i<ndim; i++) lnaxis[i] = (olong)naxis[i];
   return  ObitGPUFArrayCreate (name,(olong) ndim, (olong*)lnaxis);
}

extern void GPUFArrayToData(ObitGPUFArray* in, ObitFArray *data, ObitErr *err) {
  ObitGPUFArrayToData(in, data->array, err);
}

extern void GPUFArrayFromData(ObitGPUFArray* in, ObitFArray *data, ObitErr *err) {
  ObitGPUFArrayFromData(in, data->array, err);
}

extern void GPUFArrayToGPU(ObitGPUFArray* in, ObitErr *err) {
  ObitGPUFArrayToGPU(in, err);
}

extern void GPUFArrayToHost(ObitGPUFArray* in, ObitErr *err) {
  ObitGPUFArrayToHost(in, err);
}

 
extern ObitGPUFArray* GPUFArrayCopy  (ObitGPUFArray *in, ObitGPUFArray *out, ObitErr *err) {
  return ObitGPUFArrayCopy (in, out, err);
} // end GPUFArrayCopy 

extern ObitGPUFArray* GPUFArrayClone  (ObitGPUFArray *in, ObitGPUFArray *out, ObitErr *err) {
  ObitGPUFArrayClone (in, out, err);
  return out;
} // end GPUFArrayCopy 

extern char* GPUFArrayGetName (ObitGPUFArray* in) {
  if (!in) return "Undefined";
  return in->name;
} // end  GPUFArrayGetName

extern long GPUFArrayGetNdim (ObitGPUFArray* in) {
  return in->ndim;
} // end  GPUFArrayGetNdim

// returns an array 
extern  PyObject* GPUFArrayGetNaxis (ObitGPUFArray* in) {
  long i;
  PyObject *outList= PyList_New(in->ndim);

  for (i=0; i<in->ndim; i++) {
    PyList_SetItem(outList, i, PyInt_FromLong((long)in->naxis[i]));
  }
  return outList;
} // end  GPUFArrayGetNaxis

extern int GPUFArrayIsA (ObitFArray* in) {
  return ObitGPUFArrayIsA(in);
} // end  FGPUArrayIsA 

ObitGPUFArray* GPUFArrayRef (ObitGPUFArray* in) {
  return ObitGPUFArrayRef (in);
} // end GPUFArrayRef

ObitGPUFArray* GPUFArrayUnref (ObitGPUFArray* in) {
  if (!ObitGPUFArrayIsA(in)) return NULL;
  if (in && (in->ReferenceCount>0)) in = ObitGPUFArrayUnref (in);
  return in;
} // end GPUFArrayUnref

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitGPUFArray *me;
} GPUFArray;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitGPUFArray *me;
} GPUFArray;

%addmethods GPUFArray { 
  GPUFArray(char* name, long ndim, long *naxis) {
     GPUFArray *out;
     out = (GPUFArray *) malloc(sizeof(GPUFArray));
     if (strcmp(name, "None")) out->me = GPUFArrayCreate(name, ndim, naxis);
     else  out->me = NULL;
     return out;
   }
  ~GPUFArray() {
    self->me = GPUFArrayUnref(self->me);
    free(self);
  }
};

