/* $Id$  */  
/*--------------------------------------------------------------------*/
/* Swig module description for Image type                             */
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
#include "ObitImage.h"
#include "ObitData.h"
#include "ObitIOImageAIPS.h"
#include "ObitIOImageFITS.h"
%}


%inline %{

// size=1-> OBIT_IO_byRow, else OBIT_IO_byPlane
extern void ImageSetFITS(ObitImage *in, int size, int disk, char *file, 
                         int blc[7], int trc[7], ObitErr *err) {
  ObitIOSize fsize;

  if (size==1) fsize = OBIT_IO_byRow;
  else fsize = OBIT_IO_byPlane;
  ObitImageSetFITS(in, fsize, disk, file, blc, trc, err);
 }

extern void ImageSetAIPS(ObitImage *in, int size, int disk, int cno, int user, 
                         int blc[7], int trc[7], ObitErr *err) {
  ObitIOSize fsize;

  if (size==1) fsize = OBIT_IO_byRow;
  else fsize = OBIT_IO_byPlane;
  ObitImageSetAIPS(in, fsize, disk, cno, user, blc, trc, err);
 }

extern ObitData* ImageCastData (ObitImage* inImage) {
  return (ObitData*)inImage;
} // end ImageCastData

extern ObitImage* ImageCreate (char* name) {
  return newObitImage (name);
} // end  ImageCreate

extern ObitImage* ImageScratch (ObitImage *in, ObitErr *err) {
  return newObitImageScratch (in, err);
} // end  ImageScratch

extern PyObject* ImageInfo (ObitImage *in, ObitErr *err) {
  ObitIOImageAIPS *AIPSIO=NULL;
  ObitIOImageFITS *FITSIO=NULL;
  PyObject *outDict=PyDict_New();
  PyObject *o=NULL;

  if (err->error) return outDict;

  // Ensure in fully instantiated -assume OK if myIO exists 
  if (!in->myIO) ObitImageFullInstantiate (in, TRUE, err);
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
} // end  ImageInfo

extern ObitImage* ImageZap  (ObitImage *in, ObitErr *err) {
  return ObitImageZap (in, err);
} // end ImageZap

extern void ImageRename  (ObitImage *in, ObitErr *err) {
  ObitImageRename (in, err);
} // end ImageRename

extern ObitImage* ImageCopy  (ObitImage *in, ObitImage *out, 
			         ObitErr *err) {
  return ObitImageCopy (in, out, err);
} // end  ImageCopy

extern void ImageClone (ObitImage *in, ObitImage *out, ObitErr *err) {
   return  ObitImageClone (in, out, err);
} // end  ImageClone

extern void ImageClone2 (ObitImage *in1, ObitImage *in2, ObitImage *out, ObitErr *err) {
   return  ObitImageClone2 (in1, in2, out, err);
} // end  ImageClone

extern void ImageCloneMem (ObitImage *in, ObitImage *out, ObitErr *err) {
   return  ObitImageCloneMem (in, out, err);
} // end  ImageCloneMem

// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
// Table verion returned as outValue1
extern ObitTable* ImageGetTable (ObitImage *in, int access, 
			    char *tabType, long *outValue1, ObitErr *err) {
  ObitIOAccess laccess;
  olong loutValue1 = (olong)*outValue1;
  ObitTable *outTable=NULL;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  outTable = newObitImageTable (in, laccess, tabType, &loutValue1, err);
  *outValue1 = (long)loutValue1;
  return outTable;
} // end  ImageGetTable

extern int ImageOpen (ObitImage *in, int access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  ret = ObitImageOpen (in, laccess, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Open

// force header update 
extern void ImageDirty (ObitImage *in) {
  in->myStatus = OBIT_Modified;
} // end Dirty

extern int ImageClose (ObitImage *in, ObitErr *err) {
  ObitIOCode ret;
  ret =  ObitImageClose (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Close

extern int ImageRead (ObitImage *in, ObitErr *err) {
  ObitIOCode ret;
  in->extBuffer = FALSE;
  ret =  ObitImageRead (in, NULL, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Read

extern int ImageWrite (ObitImage *in, ObitErr *err) {
  ObitIOCode ret;
  in->extBuffer = FALSE;
  ret =  ObitImageWrite (in, NULL, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Write

extern int ImageReadFA (ObitImage *in, ObitFArray *array, 
                        ObitErr *err) {
  ObitIOCode ret;
  in->extBuffer = TRUE;
  /* Check comparability with an existing FArray */
  if (in->image) {
     if (!ObitFArrayIsCompatable(in->image, array)) {
       Obit_log_error(err, OBIT_InfoErr, "ImageReadFA: FArray incompatable");
       return 1;
    }
  }
  ret =  ObitImageRead (in, array->array, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end ReadFA

extern int ImageWriteFA (ObitImage *in, ObitFArray *array, 
                         ObitErr *err) {
  ObitIOCode ret;
  in->extBuffer = TRUE;
  /* Check comparability with an existing FArray */
  if (in->image) {
     if (!ObitFArrayIsCompatable(in->image, array)) {
       Obit_log_error(err, OBIT_InfoErr, "ImageWriteFA: FArray incompatable");
       return 1;
    }
  }
  ret =  ObitImageWrite (in, array->array, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end Write

extern int ImageGetPlane (ObitImage *in, ObitFArray *outData, int *plane,
			    ObitErr *err) {
  ObitIOCode ret;
  olong i, lplane[5];
  gboolean bad;
  ofloat *data;

  for (i=0; i<5; i++) lplane[i] = plane[i];
  if (outData) data = outData->array;  // buffer given?
  else data = NULL;
  // Check compatibility between outData and in
  bad = FALSE;
  if (outData && (in->image)) bad = !ObitFArrayIsCompatable(outData,in->image);
  else if (outData && (in->myDesc->inaxes[0]>0)) { // check dimensions if possible
    bad = (outData->ndim<2) || (outData->naxis[0]!=in->myDesc->inaxes[0]) ||
          (outData->naxis[1]!=in->myDesc->inaxes[1]);
  }
  if (bad) {
    Obit_log_error(err, OBIT_Error, 
	           "ImageGetPlane: specified FArray incompatable with image");
    return 1;
  }
  ret = ObitImageGetPlane (in, data, plane, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end ImageGetPlane

extern int ImagePutPlane (ObitImage *in, ObitFArray *outData, int *plane,
			  ObitErr *err) {
  ObitIOCode ret;
  olong i, lplane[5];
  gboolean bad;
  ofloat *data;

  for (i=0; i<5; i++) lplane[i] = plane[i];
  if (outData) data = outData->array;  // buffer given?
  else data = NULL;
  // Check compatibility between outData and in
  bad = FALSE;
  if (outData && (in->image)) bad = !ObitFArrayIsCompatable(outData,in->image);
  else if (outData && (in->myDesc->inaxes[0]>0)) { // check dimensions if possible
    bad = (outData->ndim<2) || (outData->naxis[0]!=in->myDesc->inaxes[0]) ||
          (outData->naxis[1]!=in->myDesc->inaxes[1]);
  }
  if (bad) {
    Obit_log_error(err, OBIT_Error, 
	           "ImagePutPlane: specified FArray incompatable with image");
    return 1;
  }
  ret = ObitImagePutPlane (in, data, lplane, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end ImagePutPlane

extern int ImageZapTable (ObitImage *in, char *tabType, long tabVer, 
			    ObitErr *err) {
  ObitIOCode ret;
  ret = ObitImageZapTable (in, tabType, tabVer, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  ImageZapTable

extern int ImageCopyTables (ObitImage *in, ObitImage *out, char **exclude,
		  	        char **include, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitImageCopyTables  (in, out, exclude, include, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  CopyTables

extern int ImageUpdateTables (ObitImage *in, ObitErr *err) {
  ObitIOCode ret;
  ret = ObitImageUpdateTables (in, err);
  if (ret==OBIT_IO_OK) return 0;
  else return 1;
} // end  UpdateTables

// Open and close to fully instantiate
// access 1=READONLY, 2=WRITEONLY, 3=READWRITE
extern int ImagefullInstantiate (ObitImage* in, int access, ObitErr *err) {
  ObitIOCode ret;
  ObitIOAccess laccess;

  laccess = OBIT_IO_ReadOnly;
  if (access==2) laccess = OBIT_IO_WriteOnly;
  else if (access==3) laccess = OBIT_IO_ReadWrite;
  ret = ObitImageOpen (in, laccess, err);
  ret = ObitImageClose (in, err);
  if ((err->error) || (ret!=OBIT_IO_OK)) return 1;
  else return 0;
} // end ImagefullInstantiate


extern ObitImage* ImageUnref (ObitImage* in) {
  if (!ObitImageIsA(in)) return NULL;
  if (in && (in->ReferenceCount>0)) in = ObitImageUnref(in);
  return in;
}

extern ObitImage*  ImageRef (ObitImage* in) {
  return ObitImageRef(in);
}

extern ObitInfoList* ImageGetList (ObitImage* in) {
  return ObitInfoListRef(in->info);
}

extern ObitTableList* ImageGetTableList (ObitImage* in) {
  return ObitTableListRef(in->tableList);
}

extern ObitImageDesc* ImageGetDesc (ObitImage* in) {
  return ObitImageDescRef(in->myDesc);
}

extern void ImageSetDesc (ObitImage* in, ObitImageDesc* desc) {
  in->myDesc = ObitImageDescUnref(in->myDesc);
  in->myDesc = ObitImageDescRef(desc);
}

extern ObitFArray* ImageGetFArray (ObitImage* in) {
  return ObitFArrayRef(in->image);
}

extern void ImageSetFArray (ObitImage* in, ObitFArray *image) {
  in->image = ObitFArrayUnref(in->image);
  in->image = ObitFArrayRef(image);
}

extern ObitImage* ImageGetBeam (ObitImage* in) {
  return (ObitImage*)ObitImageRef(in->myBeam);
}

extern void ImageSetBeam (ObitImage* in, ObitImage *beam) {
  in->myBeam = (Obit*)ObitImageUnref(in->myBeam);
  in->myBeam = (Obit*)ObitImageRef(beam);
}

extern long ImageGetHighVer (ObitImage* in, char *tabType) {
  return ObitTableListGetHigh(in->tableList, tabType);
}

extern int ImageisScratch (ObitImage* in) {
  return (int)in->isScratch;
}

extern int ImageIsA (ObitImage* in) {
  return ObitImageIsA(in);
}

extern char* ImageGetName (ObitImage* in) {
  if (ObitImageIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitImage *me;
} Image;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitImage *me;
} Image;

%addmethods Image { 
  Image(char* name) {
     Image *out;
     out = (Image *) malloc(sizeof(Image));
     if (strcmp(name, "None")) out->me = ImageCreate(name);
     else out->me = NULL;
     return out;
   }
  ~Image() { /* Scratch files may be deleted separately */
  if (self && (self->me)) {
    if (self->me->ReferenceCount>0) 
       self->me = ImageUnref(self->me);
    free(self);
    } 
  }
};

