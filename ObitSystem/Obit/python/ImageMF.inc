/* $Id: ImageMF.inc 36 2008-09-29 12:54:35Z bill.cotton $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageMF type                           */
/*                                                                    */
/*;  Copyright (C) 2010-2016                                          */
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


%inline %{

// size=1-> OBIT_IO_byRow, else OBIT_IO_byPlane
extern void ImageMFSetFITS(ObitImageMF *in, int size, int disk, char *file, 
                         int blc[7], int trc[7], ObitErr *err) {
  ObitIOSize fsize;

  if (size==1) fsize = OBIT_IO_byRow;
  else fsize = OBIT_IO_byPlane;
  ObitImageSetFITS(in, fsize, disk, file, blc, trc, err);
 }

extern void ImageMFSetAIPS(ObitImageMF *in, int size, int disk, int cno, int user, 
                         int blc[7], int trc[7], ObitErr *err) {
  ObitIOSize fsize;

  if (size==1) fsize = OBIT_IO_byRow;
  else fsize = OBIT_IO_byPlane;
  ObitImageSetAIPS(in, fsize, disk, cno, user, blc, trc, err);
 }

extern ObitData* ImageMFCastData (ObitImageMF* inImage) {
  return (ObitData*)inImage;
} // end ImageMFCastData

extern ObitImageMF* ImageMFCreate (char* name) {
  return newObitImageMF (name);
} // end  ImageMFCreate

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
    o = PyInt_FromLong((long)AIPSIO->disk);
    PyDict_SetItemString(outDict, "disk", o);
    o = PyInt_FromLong((long)AIPSIO->CNO);
    PyDict_SetItemString(outDict, "CNO", o);
    o = PyInt_FromLong((long)AIPSIO->UserId);
    PyDict_SetItemString(outDict, "user", o);
  } else if (ObitIOImageFITSIsA(in->myIO)) {  // FITS
    o = PyString_InternFromString("FITS");
    PyDict_SetItemString(outDict, "type", o);
    FITSIO = (ObitIOImageFITS*)in->myIO;
    o = PyInt_FromLong((long)FITSIO->disk);
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

extern ObitImageMF* ImageMFUnref (ObitImageMF* in) {
  if (!ObitImageMFIsA(in)) return NULL;
  if (in && (in->ReferenceCount>0)) in = ObitImageMFUnref(in);
  return in;
}

extern ObitImageMF*  ImageMFRef (ObitImageMF* in) {
  return ObitImageMFRef(in);
}

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

extern int ImageMFisScratch (ObitImageMF* in) {
  return (int)in->isScratch;
}

extern int ImageMFIsA (ObitImageMF* in) {
  return ObitImageMFIsA(in);
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

%addmethods ImageMF { 
  ImageMF(char* name) {
     ImageMF *out;
     out = (ImageMF *) malloc(sizeof(ImageMF));
     if (strcmp(name, "None")) out->me = ImageMFCreate(name);
     else out->me = NULL;
     return out;
   }
  ~ImageMF() { /* Scratch files may be deleted separately */
  if (self && (self->me)) {
    if (self->me->ReferenceCount>0) 
       self->me = ImageMFUnref(self->me);
    free(self);
    } 
  }
};

