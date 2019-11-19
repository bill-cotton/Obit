/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for PixHistFDR type                        */
/*                                                                    */
/*;  Copyright (C) 2011,2019                                          */
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
#include "ObitPixHisto.h"
#include "ObitImage.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitPixHisto *me;
} PixHistFDR;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitPixHisto *me;
} PixHistFDR;


%inline %{
extern ObitPixHisto* newPixHistFDR (char* name) {
  return newObitPixHisto (name);
} // end  newPixHistFDR

extern ObitPixHisto* PixHistFDRUnref (ObitPixHisto* in) {
  if (!ObitPixHistoIsA(in)) return NULL;
  return ObitPixHistoUnref(in);
}

extern ObitPixHisto*  PixHistFDRRef (ObitPixHisto* in) {
  return ObitPixHistoRef(in);
}

extern ObitPixHisto* PixHistFDRCreate (char *name, ObitImage *image, ObitErr *err) {
  return ObitPixHistoCreate((gchar*)name, image, err);
}

void CreateOPixHistFDR(PyObject *self, char* name) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_PixHistFDR, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Create: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((PixHistFDR*)ptr)->me = newPixHistFDR(name);
   else ((PixHistFDR*)ptr)->me = NULL;
   return;
 } // end CreateOPixHistFDR

void DeleteOPixHistFDR(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_PixHistFDR, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Delete: could not recover c struct");
	return;
   }
   ((PixHistFDR*)ptr)->me = PixHistFDRUnref(((PixHistFDR*)ptr)->me);
} // end DeleteOPixHistFDR

ObitPixHisto* OPixHistFDR_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_PixHistFDR, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"Get_me: could not recover c struct");
	return NULL;
   }
   return ((PixHistFDR*)ptr)->me;
} // end OPixHistFDR_Get_me

void OPixHistFDR_Set_me(PyObject *self, ObitPixHisto *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_PixHistFDR, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"PixHistFDR_Set_me: could not recover c struct");
	return;
   }
   ((PixHistFDR*)ptr)->me =  PixHistFDRRef(me);
}  // end OPixHistFDR_Set_me

extern ObitPixHisto* PixHistFDRCopy  (ObitPixHisto *in, ObitPixHisto *out, 
				    ObitErr *err) {
  return ObitPixHistoCopy (in, out, err);
} // end  PixHistFDRCopy

// Image
extern ObitFArray* PixHistFDRGetImage (ObitPixHisto* in) {
  return ObitImageRef(in->image);
}

// Image pixel array
extern ObitFArray* PixHistFDRGetImPix (ObitPixHisto* in) {
  return ObitFArrayRef(in->imagePix);
}

extern ObitInfoList* PixHistFDRGetList (ObitPixHisto* in) {
  return ObitInfoListRef(in->info);
}

extern void PixHistFDRHisto (ObitPixHisto* in, long *blc, long *trc, 
                             long nbin, float range, ObitErr *err) {
  olong lblc[2], ltrc[2];

  lblc[0] = (olong)blc[0]; lblc[1] = (olong)blc[1]; 
  ltrc[0] = (olong)trc[0]; ltrc[1] = (olong)trc[1]; 
  ObitPixHistoHisto(in, lblc, ltrc, (olong)nbin, (ofloat)range, err);
}

// Differential histogram
extern ObitFArray* PixHistFDRGetHisto (ObitPixHisto* in) {
  return ObitFArrayRef(in->histo);
}

// Integrated histogram
extern ObitFArray* PixHistFDRGetIntHisto (ObitPixHisto* in) {
  return ObitFArrayRef(in->intHisto);
}

extern float PixHistFDRSigma (ObitPixHisto* in) {
  return (float)in->sigma;
}


extern float PixHistFDRFlux (ObitPixHisto* in, float maxFDR, 
	                  ObitErr *err) {
 return (float)ObitPixHistoFDRFlux(in, (ofloat)maxFDR, err);
}

extern char* PixHistFDRGetName (ObitPixHisto* in) {
  if (ObitPixHistoIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern long PixHistFDRIsA (ObitPixHisto* in) {
  return (long)ObitPixHistoIsA(in);
}

%}

