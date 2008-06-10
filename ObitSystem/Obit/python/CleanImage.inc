/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for CleanImage type                        */
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
#include "ObitDConCleanImage.h"
%}


%inline %{
extern ObitDConCleanImage* newCleanImage (char* name) {
  return newObitDConCleanImage (name);
} // end  newCleanImage

extern ObitDConCleanImage* CleanImageCopy  (ObitDConCleanImage *in, 
	ObitDConCleanImage *out, 
        ObitErr *err) {
  return ObitDConCleanImageCopy (in, out, err);
} // end  CleanImageCopy

extern ObitDConCleanImage* CleanImageUnref (ObitDConCleanImage* in) {
  if (!ObitDConCleanImageIsA(in)) return NULL;
  return ObitDConCleanImageUnref(in);
}

extern ObitDConCleanImage*  CleanImageRef (ObitDConCleanImage* in) {
  return ObitDConCleanImageRef(in);
}

extern ObitInfoList* CleanImageGetList (ObitDConCleanImage* in) {
  return ObitInfoListRef(in->info);
}

extern ObitImageMosaic* CleanImageGetImageMosaic (ObitDConCleanImage* in) {
  return ObitImageMosaicRef(in->mosaic);
}

extern void CleanImageSetImageMosaic (ObitDConCleanImage* in, ObitImageMosaic *mosaic) {
  in->mosaic = ObitImageMosaicUnref(in->mosaic);  /* Out with the old */
  in->mosaic = ObitImageMosaicRef(mosaic);        /* In with the new */
}

extern ObitDConCleanWindow* CleanImageGetWindow (ObitDConCleanImage* in) {
  return ObitDConCleanWindowRef(in->window);
}

extern void CleanImageSetWindow (ObitDConCleanImage* in, ObitDConCleanWindow *window) {
  in->mosaic = ObitDConCleanWindowUnref(in->window);  /* Out with the old */
  in->mosaic = ObitDConCleanWindowRef(window);        /* In with the new */
}

// if (win[0]<0 then this is a round window then win[1]=radius, [2,3] = center
// else rectangular and blc=(win[0],win[1]), trc= blc=(win[2],win[3])
extern void CleanImageAddWindow (ObitDConCleanImage* in, int field, int *win, 
                                 ObitErr *err) {
  olong lfield, window[4];
  ObitDConCleanWindowType type;

  // which type? 
  if (win[0]<0) { // circle
    type = OBIT_DConCleanWindow_round;
    window[0] = win[1];
    window[1] = win[2];
    window[2] = win[3];
  } else { // rectangular
    type = OBIT_DConCleanWindow_rectangle;
    window[0] = win[0];
    window[1] = win[1];
    window[2] = win[2];
    window[3] = win[3];
  }
	
  lfield = field;
  ObitDConCleanWindowAdd (in->window, lfield, type, window, err);
}  // end CleanImageAddWindow 

extern ObitDConCleanImage* CleanImageCreate (char *name, ObitImageMosaic* mosaic, 
                                             ObitErr *err) {
 return ObitDConCleanImageCreate(name, mosaic, err);
}

extern void CleanImageDeconvolve (ObitDConCleanImage* in, ObitErr *err) {
 ObitDConCleanImageDeconvolve((ObitDCon*)in, err);
}

extern void CleanImageDefWindow (ObitDConCleanImage* in, ObitErr *err) {
 ObitDConCleanDefWindow((ObitDConClean*)in, err);
}

extern char* CleanImageGetName (ObitDConCleanImage* in) {
  if (ObitDConCleanImageIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int CleanImageIsA (ObitDConCleanImage* in) {
  return ObitDConCleanImageIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitDConCleanImage *me;
} CleanImage;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitDConCleanImage *me;
} CleanImage;

%addmethods CleanImage { 
  CleanImage(char* name) {
     CleanImage *out;
     out = (CleanImage *) malloc(sizeof(CleanImage));
     if (strcmp(name, "None")) out->me = newCleanImage(name);
     else out->me = NULL;
     return out;
   }
  ~CleanImage() {
   if (self->me-> ReferenceCount>0) 
      self->me = CleanImageUnref(self->me);
   free(self);
  }
};

