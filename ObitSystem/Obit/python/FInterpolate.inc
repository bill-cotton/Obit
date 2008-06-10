/* $Id: FInterpolate.inc,v 1.3 2005/02/06 02:00:38 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for FInterpolate type                      */
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
#include "ObitFInterpolate.h"
%}

%inline %{
extern ObitFInterpolate* 
FInterpolateCreate (char* name, ObitFArray *array, ObitImageDesc *desc, long hwidth) {
  return newObitFInterpolateCreate (name, array, desc, hwidth);
} // end FInterpolateCreate

extern ObitFInterpolate* FInterpolateCopy  (ObitFInterpolate *in, ObitFInterpolate *out, 
			   ObitErr *err) {
  return ObitFInterpolateCopy (in, out, err);
} // end FInterpolateCopy


extern ObitFInterpolate* FInterpolateClone (ObitFInterpolate *in, ObitFInterpolate *out) {
  return ObitFInterpolateClone (in, out);
} // end FInterpolateClone 

extern void FInterpolateReplace (ObitFInterpolate *in, ObitFArray *newArray) {
  ObitFInterpolateReplace (in, newArray);
} // end FInterpolateReplac

extern float FInterpolatePixel (ObitFInterpolate *in, float *pixel, ObitErr *err) {
  return ObitFInterpolatePixel (in, pixel, err);
} // end FInterpolatePixel 

extern float FInterpolate1D (ObitFInterpolate *in, float pixel) {
  return ObitFInterpolate1D (in, pixel);
} // end FInterpolate1D

extern float FInterpolatePosition (ObitFInterpolate *in, double *coord, ObitErr *err) {
  return ObitFInterpolatePosition (in, coord, err);
} // end FInterpolatePosition

extern ObitInfoList* FInterpolateGetList (ObitFInterpolate* in) {
  return ObitInfoListRef(in->info);
}

extern ObitFArray* FInterpolateGetFArray (ObitFInterpolate* in) {
  return ObitFArrayRef(in->myArray);
}

extern ObitImageDesc* FInterpolateGetDesc (ObitFInterpolate* in) {
  return ObitImageDescRef(in->myDesc);
}

extern void FInterpolateSetDesc (ObitFInterpolate* in, ObitImageDesc *desc) {
  in->myDesc = ObitImageDescUnref(in->myDesc);
  in->myDesc = ObitImageDescRef(desc);
}

extern long FInterpolateGetHwidth (ObitFInterpolate* in) {
  return in->hwidth;
}

extern void FInterpolateSetHwidth (ObitFInterpolate* in, long hwidth) {
  in->hwidth = hwidth;
}

ObitFInterpolate* FInterpolateRef (ObitFInterpolate* in) {
  return ObitFInterpolateRef (in);
} // end FInterpolateRef

ObitFInterpolate* FInterpolateUnref (ObitFInterpolate* in) {
  if (!ObitFInterpolateIsA(in)) return NULL;
  return ObitFInterpolateUnref (in);
} // end FInterpolateUnref

extern int FInterpolateIsA (ObitFInterpolate* in) {
  return ObitFInterpolateIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitFInterpolate *me;
} FInterpolate;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitFInterpolate *me;
} FInterpolate;

%addmethods FInterpolate { 
  FInterpolate(char* name, ObitFArray *array, ObitImageDesc *desc, long hwidth) {
     FInterpolate *out;
     /* just create  Python structure here */
     out = (FInterpolate *) malloc(sizeof(FInterpolate));
     if (strcmp(name, "None")) out->me = FInterpolateCreate(name, array, desc, hwidth);
     else out->me = NULL;
     return out;
   }
  ~FInterpolate() {
    self->me = FInterpolateUnref(self->me);
    free(self);
  }
};


