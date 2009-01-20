/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageFit type                          */
/*                                                                    */
/*;  Copyright (C) 2007                                               */
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
#include "ObitImageFit.h"
#include "ObitImage.h"
#include "ObitFitRegion.h"
%}


%inline %{
extern ObitImageFit* newImageFit (char* name) {
  return newObitImageFit (name);
} // end  newImageFit

extern ObitImageFit* ImageFitCopy  (ObitImageFit *in, ObitImageFit *out, 
				    ObitErr *err) {
  return ObitImageFitCopy (in, out, err);
} // end  ImageFitCopy

extern ObitImageFit* ImageFitUnref (ObitImageFit* in) {
  if (!ObitImageFitIsA(in)) return NULL;
  return ObitImageFitUnref(in);
}

extern ObitImageFit*  ImageFitRef (ObitImageFit* in) {
  return ObitImageFitRef(in);
}

extern ObitInfoList* ImageFitGetList (ObitImageFit* in) {
  return ObitInfoListRef(in->info);
}

extern int ImageFitFit (ObitImageFit* in, ObitImage *image, 
                        ObitFitRegion* reg, ObitErr *err) {
 return (int)ObitImageFitFit(in, image, reg, err);
}

extern char* ImageFitGetName (ObitImageFit* in) {
  if (ObitImageFitIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int ImageFitIsA (ObitImageFit* in) {
  return ObitImageFitIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitImageFit *me;
} ImageFit;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitImageFit *me;
} ImageFit;

%addmethods ImageFit { 
  ImageFit(char* name) {
     ImageFit *out;
     out = (ImageFit *) malloc(sizeof(ImageFit));
     if (strcmp(name, "None")) {
        out->me = newImageFit(name);
     } else out->me = NULL;
     return out;
   }
  ~ImageFit() {
   if (!self) return;  // Not defined
   if (self->me->ReferenceCount>0) 
      self->me = ImageFitUnref(self->me);
   free(self);
  }
};

