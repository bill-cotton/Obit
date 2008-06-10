/* $Id$   */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageMosaic type                       */
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
#include "ObitImageMosaic.h"
%}


%inline %{
extern ObitImageMosaic* newImageMosaic (char* name, int number) {
  return newObitImageMosaic (name, number);
} // end  newImageMosaic

extern ObitImageMosaic* ImageMosaicCopy  (ObitImageMosaic *in, ObitImageMosaic *out, 
				          ObitErr *err) {
  return ObitImageMosaicCopy (in, out, err);
} // end  ImageMosaicCopy

extern ObitImageMosaic* ImageMosaicUnref (ObitImageMosaic* in) {
  if (!ObitImageMosaicIsA(in)) return NULL;
  return ObitImageMosaicUnref(in);
}

extern ObitImageMosaic*  ImageMosaicRef (ObitImageMosaic* in) {
  return ObitImageMosaicRef(in);
}

extern void  ImageMosaicZapImage(ObitImageMosaic* in, int number, 
                                 ObitErr *err ) {
  return ObitImageMosaicZapImage(in, number, err);
}

extern ObitInfoList* ImageMosaicGetList (ObitImageMosaic* in) {
  return ObitInfoListRef(in->info);
}

extern ObitImage* ImageMosaicGetFullImage (ObitImageMosaic* in, ObitErr *err) {
  return ObitImageMosaicGetFullImage(in, err);
}

extern void ImageMosaicSetFullImage (ObitImageMosaic* in, 
                                 ObitImage *image, ObitErr *err) {
 ObitImageMosaicSetFullImage(in, image, err);
}

extern ObitImage* ImageMosaicGetImage (ObitImageMosaic* in, int number, 
                                       ObitErr *err) {
  return ObitImageMosaicGetImage(in, number, err);
}

extern void ImageMosaicSetImage (ObitImageMosaic* in, int number, 
                                 ObitImage *image, ObitErr *err) {
 ObitImageMosaicSetImage(in, number, image, err);
}

extern ObitImageMosaic* ImageMosaicCreate (char *name, ObitUV *uvData, ObitErr *err) {
 return ObitImageMosaicCreate(name, uvData, err);
}

extern void ImageMosaicDefine (ObitImageMosaic* in, ObitUV *uvData, int doBeam,
                               ObitErr *err) {
 gboolean LdoBeam = doBeam;
 ObitImageMosaicDefine(in, uvData, LdoBeam, err);
}

extern void ImageMosaicFlatten (ObitImageMosaic* in, ObitErr *err) {
 ObitImageMosaicFlatten(in, err);
}

extern char* ImageMosaicGetName (ObitImageMosaic* in) {
  if (ObitImageMosaicIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int ImageMosaicGetNumber (ObitImageMosaic* in) {
  if (ObitImageMosaicIsA(in)) {
    return in->numberImages;
  } else {
    return 0;
  }
}

extern int ImageMosaicIsA (ObitImageMosaic* in) {
  return ObitImageMosaicIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitImageMosaic *me;
} ImageMosaic;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitImageMosaic *me;
} ImageMosaic;

%addmethods ImageMosaic { 
  ImageMosaic(char* name, int number) {
     ImageMosaic *out;
     out = (ImageMosaic *) malloc(sizeof(ImageMosaic));
     if (strcmp(name, "None")) out->me = newImageMosaic(name, number);
     else out->me = NULL;
     return out;
   }
  ~ImageMosaic() {
   if (self->me-> ReferenceCount>0) 
      self->me = ImageMosaicUnref(self->me);
   free(self);
  }
};

