/* $Id:  $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageMF type                           */
/*                                                                    */
/*;  Copyright (C) 2010-2019                                          */
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
#include "ObitImageMF.h"
#include "ObitData.h"
#include "ObitIOImageAIPS.h"
#include "ObitIOImageFITS.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitImageMF *me;
} ImageMF;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitImageMF *me;
} ImageMF;

%inline %{
extern ObitImageMF* ImageMFCreate (char* name) {
  return newObitImageMF (name);
} // end  ImageMFCreate

extern ObitImageMF* ImageMFUnref (ObitImageMF* in) {
  if (!ObitImageMFIsA(in)) return NULL;
  if (in && (in->ReferenceCount>0)) in = ObitImageMFUnref(in);
  return in;
}

extern ObitImageMF*  ImageMFRef (ObitImageMF* in) {
  return ObitImageMFRef(in);
}

void CreateImageMF(PyObject *self, char* name) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ImageMF, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateImageMF: could not recover c struct");
	return;
   }
   if (strcmp(name, "None")) ((ImageMF*)ptr)->me = ImageMFCreate(name);
   else  ((ImageMF*)ptr)->me = NULL;
   return;
 } // end CreateImageMF

void DeleteImageMF(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ImageMF, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteImageMF: could not recover c struct");
	return;
   }
   ((ImageMF*)ptr)->me = ObitImageMFUnref(((ImageMF*)ptr)->me);
}// end DeleteImageMF

ObitImageMF* ImageMF_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ImageMF, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"ImageMF_Get_me: could not recover c struct");
	return NULL;
   }
   return ((ImageMF*)ptr)->me;
}  // end ImageMF_Get_me

void ImageMF_Set_me(PyObject *self, ObitImageMF *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_ImageMF, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"ImageMF_Set_me: could not recover c struct");
	return;
   }
   ((ImageMF*)ptr)->me =  ObitImageMFRef(me);
}  // end ImageMF_Set_me

// size=1-> OBIT_IO_byRow, else OBIT_IO_byPlane
extern void ImageMFSetFITS(ObitImageMF *in, long size, long disk, char *file, 
                         long *blc, long *trc, ObitErr *err) {
  ObitIOSize fsize;
  olong i, lblc[7], ltrc[7], ldisk;
  if (size==1) fsize = OBIT_IO_byRow;
  else fsize = OBIT_IO_byPlane;
  for (i=0; i<7; i++) {lblc[i]=(olong)blc[i]; ltrc[i]=(olong)trc[i];}
  ldisk = (olong)disk;
  ObitImageSetFITS(in, fsize, ldisk, file, lblc, ltrc, err);
 }

extern void ImageMFSetAIPS(ObitImageMF *in, long size, long disk, long cno, long user, 
                         long *blc, long *trc, ObitErr *err) {
  ObitIOSize fsize;
  olong i, lblc[7], ltrc[7], ldisk, lcno, luser;
  if (size==1) fsize = OBIT_IO_byRow;
  else fsize = OBIT_IO_byPlane;
  for (i=0; i<7; i++) {lblc[i]=(olong)blc[i]; ltrc[i]=(olong)trc[i];}
  ldisk = (olong)disk; lcno = (olong)cno; luser = (olong)user;
  ObitImageSetAIPS(in, fsize, ldisk, lcno, luser, lblc, ltrc, err);
 }

extern ObitData* ImageMFCastData (ObitImageMF* inImage) {
  return (ObitData*)inImage;
} // end ImageMFCastData

extern ObitImage* ImageMFCastImage (ObitImageMF* inImage) {
  return (ObitImage*)inImage;
} // end ImageMFCastImage

// Open and close to fully instantiate
// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
extern long ImageMFfullInstantiate (ObitImageMF* in, long access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  ret = ObitImageOpen ((ObitImage*)in, laccess, err);
  ret = ObitImageClose ((ObitImage*)in, err);
  if ((err->error) || (ret!=OBIT_IO_OK)) return 1;
  else return 0;
} // end ImageMFfullInstantiate

extern PyObject* ImageMFInfo (ObitImageMF *in, ObitErr *err) {
  ObitIOImageAIPS *AIPSIO=NULL;
  ObitIOImageFITS *FITSIO=NULL;
  PyObject *outDict=PyDict_New();
  PyObject *o=NULL;

  if (err->error) return outDict;

  // Ensure in fully instantiated -assume OK if myIO exists 
  if (!in->myIO) ObitImageFullInstantiate ((ObitImage*)in, TRUE, err);
  if (err->error) return outDict;

  // Get details and save in dict
  if (ObitIOImageAIPSIsA(in->myIO)) {  // AIPS
    o = PyString_InternFromString("AIPS");
    PyDict_SetItemString(outDict, "type", o);
    AIPSIO = (ObitIOImageAIPS*)in->myIO;
    o = PyLong_FromLong((long)AIPSIO->disk);
    PyDict_SetItemString(outDict, "disk", o);
    o = PyLong_FromLong((long)AIPSIO->CNO);
    PyDict_SetItemString(outDict, "CNO", o);
    o = PyLong_FromLong((long)AIPSIO->UserId);
    PyDict_SetItemString(outDict, "user", o);
  } else if (ObitIOImageFITSIsA(in->myIO)) {  // FITS
    o = PyString_InternFromString("FITS");
    PyDict_SetItemString(outDict, "type", o);
    FITSIO = (ObitIOImageFITS*)in->myIO;
    o = PyLong_FromLong((long)FITSIO->disk);
    PyDict_SetItemString(outDict, "disk", o);
    o = PyString_InternFromString((char*)FITSIO->FileName);
    PyDict_SetItemString(outDict, "filename", o);
  } else {  // Don't know this one
    o = PyString_InternFromString("UNKNOWN");
    PyDict_SetItemString(outDict, "type", o);
  }
  return outDict;
} // end  ImageMFInfo

extern ObitImageMF* ImageMFZap  (ObitImageMF *in, ObitErr *err) {
  ObitImageMFZap ((ObitImage*)in, err);
  return NULL;
} // end ImageMFZap

extern void ImageMFRename  (ObitImageMF *in, ObitErr *err) {
  ObitImageRename ((ObitImage*)in, err);
} // end ImageMFRename

extern ObitImageMF* ImageMFCopy  (ObitImageMF *in, ObitImageMF *out, 
			         ObitErr *err) {
  return ObitImageMFCopy (in, out, err);
} // end  ImageMFCopy

// force header update 
extern void ImageMFDirty (ObitImageMF *in) {
  in->myStatus = OBIT_Modified;
} // end Dirty

extern ObitInfoList* ImageMFGetList (ObitImageMF* in) {
  return ObitInfoListRef(in->info);
}

extern ObitTableList* ImageMFGetTableList (ObitImageMF* in) {
  return ObitTableListRef(in->tableList);
}

extern ObitImageDesc* ImageMFGetDesc (ObitImageMF* in) {
  return ObitImageDescRef(in->myDesc);
}

extern void ImageMFSetDesc (ObitImageMF* in, ObitImageDesc* desc) {
  in->myDesc = ObitImageDescUnref(in->myDesc);
  in->myDesc = ObitImageDescRef(desc);
}

extern ObitFArray* ImageMFGetFArray (ObitImageMF* in) {
  return ObitFArrayRef(in->image);
}

extern void ImageMFSetFArray (ObitImageMF* in, ObitFArray *image) {
  in->image = ObitFArrayUnref(in->image);
  in->image = ObitFArrayRef(image);
}

extern ObitImageMF* ImageMFGetBeam (ObitImageMF* in) {
  return (ObitImageMF*)ObitImageMFRef(in->myBeam);
}

extern void ImageMFSetBeam (ObitImageMF* in, ObitImageMF *beam) {
  in->myBeam = (Obit*)ObitImageMFUnref(in->myBeam);
  in->myBeam = (Obit*)ObitImageMFRef(beam);
}

extern long ImageMFGetHighVer (ObitImageMF* in, char *tabType) {
  return ObitTableListGetHigh(in->tableList, tabType);
}

extern long ImageMFisScratch (ObitImageMF* in) {
  return (long)in->isScratch;
}

extern long ImageMFIsA (ObitImageMF* in) {
  return (long)ObitImageMFIsA(in);
}

extern char* ImageMFGetName (ObitImageMF* in) {
  if (ObitImageMFIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern void ImageMFFitSpec (ObitImageMF* in, float antSize, ObitErr* err) {
  ObitImageMFFitSpec (in, antSize, err);
}

extern void ImageMFFitSpec2 (ObitImageMF* in, ObitImageMF* out, ObitErr* err) {
  ObitImageMFFitSpec2 (in, out, err);
}

void ObitImageMFEffFqCorr (ObitImageMF *in, ObitImageMF *in2, ObitErr *err);
extern void ImageMFEffFqCorr (ObitImageMF* in, ObitImageMF* in2, ObitErr* err) {
  ObitImageMFEffFqCorr (in, in2, err);
}

%}

