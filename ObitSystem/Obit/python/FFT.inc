/* $Id$  */  
/*--------------------------------------------------------------------*/
/* Swig module description for ObitFFT type                           */
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
#include "ObitFFT.h"
%}


%inline %{
// dir   1 => OBIT_FFT_Forward(R2C), 2 => OBIT_FFT_Reverse(C2R)
// type  1 => OBIT_FFT_FullComplex,  2 => OBIT_FFT_HalfComplex
extern ObitFFT* FFTCreate(char* name, int dir, int type, int rank, int *dim) {
   ObitFFTdir ldir;
   ObitFFTtype ltype;
   if (dir==1) ldir = OBIT_FFT_Forward;
   else ldir = OBIT_FFT_Reverse;
   if (type==1) ltype = OBIT_FFT_FullComplex;
   else ltype = OBIT_FFT_HalfComplex;
   return newObitFFT  (name, ldir, ltype, rank, dim);
}

extern int FFTSuggestSize  (int length) {
  return ObitFFTSuggestSize(length);
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

extern int FFTGetRank (ObitFFT* in) {
  return in->rank;
} // end  FFTGetRank

// returns an array of 7 elements no matter
extern int* FFTGetDim (ObitFFT* in) {
  return in->dim;
} // end  FFTGetDim

extern int FFTIsA (ObitFFT* in) {
  return ObitFFTIsA(in);
} // end  FFTIsA 

ObitFFT* FFTRef (ObitFFT* in) {
  return ObitFFTRef (in);
} // end FFTRef

ObitFFT* FFTUnref (ObitFFT* in) {
  if (!ObitFFTIsA(in)) return NULL;
  return ObitFFTUnref (in);
} // end FFTUnref



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

%addmethods FFT { 
  FFT(char* name, int dir, int type, int rank, int *dim) {
     FFT *out;
     /* just create  Python structure here */
     out = (FFT *) malloc(sizeof(FFT));
     if (strcmp(name, "None")) out->me = FFTCreate(name, dir, type, rank, dim);
     else out->me = NULL;
     return out;
   }
  ~FFT() {
    self->me = FFTUnref(self->me);
    free(self);
  }
};

