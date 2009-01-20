/* $Id$ */  
/*--------------------------------------------------------------------*/
/* Swig module description for CleanOTF type                          */
/*                                                                    */
/*;  Copyright (C) 2006                                               */
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
#include "ObitDConCleanOTFRec.h"
%}

%inline %{
extern ObitDConCleanOTFRec* CleanOTFRecCreate (char* name, 
			 ObitOTF *inOTF, ObitErr *err) {
  return ObitDConCleanOTFRecCreate (name, inOTF, err);
} // end  CleanOTFRecCreate

extern void CleanOTFRecClean (ObitDConCleanOTFRec *in, ObitErr *err) {
  ObitDConCleanOTFRecDeconvolve ((ObitDCon*)in, err);
} // end CleanOTFRecClean

extern void CleanOTFRecRestore (ObitDConCleanOTFRec *in, ObitErr *err) {
  ObitDConCleanOTFRecRestore ((ObitDConClean*)in, err);
} // end CleanOTFRecRestore

extern ObitDConCleanOTFRec* CleanOTFRecUnref (ObitDConCleanOTFRec* in) {
  if (!ObitDConCleanOTFRecIsA(in)) return NULL;
  return ObitDConCleanOTFRecUnref(in);
}

extern ObitDConCleanOTFRec*  CleanOTFRecRef (ObitDConCleanOTFRec* in) {
  return ObitDConCleanOTFRecRef(in);
}

extern ObitImageMosaic* CleanOTFRecGetMosaic (ObitDConCleanOTFRec* in) {
  return ObitImageMosaicRef(in->mosaic);
}

extern ObitOTF* CleanOTFRecGetOTF (ObitDConCleanOTFRec* in) {
  return ObitOTFRef(in->myOTF);
}

extern ObitImage* CleanOTFRecGetBeam (ObitDConCleanOTFRec* in) {
  return ObitImageRef(in->beam);
}

extern void CleanOTFRecSetBeam (ObitDConCleanOTFRec* in, ObitImage* image) {
  in->beam = ObitImageUnref(in->beam);
  in->beam = ObitImageRef(image);
}

extern ObitImage* CleanOTFRecGetClean (ObitDConCleanOTFRec* in) {
  return ObitImageRef(in->clean);
}

extern void CleanOTFRecSetClean (ObitDConCleanOTFRec* in, ObitImage* image) {
  in->clean = ObitImageUnref(in->clean);
  in->clean = ObitImageRef(image);
}

extern ObitImage* CleanOTFRecGetWeight (ObitDConCleanOTFRec* in) {
  return ObitImageRef(in->weight);
}

extern void CleanOTFRecSetWeight (ObitDConCleanOTFRec* in, ObitImage* image) {
  in->weight = ObitImageUnref(in->weight);
  in->weight = ObitImageRef(image);
}

extern void CleanOTFRecDefWindow (ObitDConCleanOTFRec* in, ObitErr *err) {
 ObitDConCleanDefWindow((ObitDConClean*)in, err);
}

extern ObitDConCleanWindow* CleanOTFRecGetWindow (ObitDConCleanOTFRec* in) {
  return ObitDConCleanWindowRef(in->window);
}

extern void CleanOTFRecSetWindow (ObitDConCleanOTFRec* in, ObitDConCleanWindow *window) {
  in->mosaic = ObitDConCleanWindowUnref(in->window);  /* Out with the old */
  in->mosaic = ObitDConCleanWindowRef(window);        /* In with the new */
}

// if (win[0]<0 then this is a round window then win[1]=radius, [2,3] = center
// else rectangular and blc=(win[0],win[1]), trc= blc=(win[2],win[3])
extern void CleanOTFRecAddWindow (ObitDConCleanOTFRec* in, int *win, 
                                 ObitErr *err) {
  olong window[4];
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
	
  ObitDConCleanWindowAdd (in->window, 1, type, window, err);
}  // end CleanOTFRecAddWindow 

extern ObitInfoList* CleanOTFRecGetList (ObitDConCleanOTFRec* in) {
  return ObitInfoListRef(in->info);
}

extern long CleanOTFRecGetNiter (ObitDConCleanOTFRec* in) {
  return in->niter;
}

extern float CleanOTFRecGetGain (ObitDConCleanOTFRec* in) {
  return in->gain[0];
}

extern float CleanOTFRecGetFlux (ObitDConCleanOTFRec* in) {
  return in->minFlux[0];
}

extern float CleanOTFRecGetCleanSize (ObitDConCleanOTFRec* in) {
  return in->cleanSize;
}

extern int CleanOTFRecIsA (ObitDConCleanOTFRec* in) {
  return ObitDConCleanOTFRecIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitDConCleanOTFRec *me;
} CleanOTFRec;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitDConCleanOTFRec *me;
} CleanOTFRec;

%addmethods CleanOTFRec { 
  CleanOTFRec(char* name, ObitOTF *inOTF, ObitErr *err) {
     CleanOTFRec *out;
     out = (CleanOTFRec *) malloc(sizeof(CleanOTFRec));
     if(name) out->me = CleanOTFRecCreate(name, inOTF, err);
     else  out->me = NULL;
     return out;
   }
  ~CleanOTFRec() {
    self->me = CleanOTFRecUnref(self->me);
    free(self);
  }
};



