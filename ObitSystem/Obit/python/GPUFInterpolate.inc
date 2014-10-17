/* $Id:  $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for GPUFInterpolate type                   */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

%{
#include "ObitGPUFInterpolate.h"
%}

%inline %{
extern ObitGPUFInterpolate* 
GPUFInterpolateCreate (char* name, ObitFArray *inArray,  
                       ObitFArray *xArray, ObitFArray *yArray, 
		       long hwidth, ObitErr *err) {
  return ObitGPUFInterpolateCreate (name, inArray, xArray, yArray, hwidth, err);
} // end GPUFInterpolateCreate

extern void 
GPUFInterpolateImage (ObitGPUFInterpolate *in, ObitFArray *inArray,
		      ObitFArray *outArray, ObitErr *err) {
  ObitGPUFInterpolateImage (in, inArray, outArray, err);
} // end GPUFInterpolateImage

extern ObitGPUFInterpolate* 
GPUFInterpolateCopy  (ObitGPUFInterpolate *in, ObitGPUFInterpolate *out, 
		      ObitErr *err) {
  return ObitGPUFInterpolateCopy (in, out, err);
} // end GPUFInterpolateCopy

extern void 
GPUFInterpolateClone (ObitGPUFInterpolate *in, ObitGPUFInterpolate *out, ObitErr *err) {
  ObitGPUFInterpolateClone (in, out, err);
} // end GPUFInterpolateClone 

ObitGPUFInterpolate* GPUFInterpolateRef (ObitGPUFInterpolate* in) {
  return ObitGPUFInterpolateRef (in);
} // end GPUFInterpolateRef

ObitGPUFInterpolate* GPUFInterpolateUnref (ObitGPUFInterpolate* in) {
  if (!ObitGPUFInterpolateIsA(in)) return NULL;
  return ObitGPUFInterpolateUnref (in);
} // end GPUFInterpolateUnref

extern int GPUFInterpolateIsA (ObitGPUFInterpolate* in) {
  return ObitGPUFInterpolateIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitGPUFInterpolate *me;
} GPUFInterpolate;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitGPUFInterpolate *me;
} GPUFInterpolate;

%addmethods GPUFInterpolate { 
  GPUFInterpolate(char* name, ObitFArray *inArray, 
                  ObitFArray *xArray, ObitFArray *yArray, long hwidth, ObitErr *err) {
     GPUFInterpolate *out;
     /* just create  Python structure here */
     out = (GPUFInterpolate *) malloc(sizeof(GPUFInterpolate));
     if (strcmp(name, "None")) 
        out->me = GPUFInterpolateCreate(name, inArray, xArray, yArray, hwidth, err);
     else out->me = NULL;
     return out;
   }
  ~GPUFInterpolate() {
    self->me = GPUFInterpolateUnref(self->me);
    free(self);
  }
};


