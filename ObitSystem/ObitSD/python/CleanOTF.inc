/* $Id: CleanOTF.inc,v 1.4 2005/08/07 18:48:56 bcotton Exp $ */  
/*--------------------------------------------------------------------*/
/* Swig module description for CleanOTF type                          */
/*                                                                    */
/*;  Copyright (C) 2005                                               */
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
#include "ObitDConCleanOTF.h"
%}

%inline %{
extern ObitDConCleanOTF* CleanOTFCreate (char* name, 
			 ObitImage *dirty, ObitImage *beam,
			 ObitImage *clean, ObitErr *err) {
// May not be passed a beam, name of beam will be "NoBeam"
   ObitImage *lbeam=NULL;
   if (beam && strncmp(beam->name,"NoBeam",6)) lbeam = beam;
  return ObitDConCleanOTFCreate (name, dirty, lbeam, clean, err);
} // end  CleanOTFCreate

extern void CleanOTFClean (ObitDConCleanOTF *in, ObitErr *err) {
  ObitDConCleanOTFDeconvolve ((ObitDCon*)in, err);
} // end CleanOTFClean

extern void CleanOTFRestore (ObitDConCleanOTF *in, ObitErr *err) {
  ObitDConCleanOTFRestore ((ObitDConClean*)in, err);
} // end CleanOTFRestore

extern ObitDConCleanOTF* CleanOTFUnref (ObitDConCleanOTF* in) {
  if (!ObitDConCleanOTFIsA(in)) return NULL;
  return ObitDConCleanOTFUnref(in);
}

extern ObitDConCleanOTF*  CleanOTFRef (ObitDConCleanOTF* in) {
  return ObitDConCleanOTFRef(in);
}

extern ObitImage* CleanOTFGetDirty (ObitDConCleanOTF* in) {
  return ObitImageRef(in->dirty);
}

extern void CleanOTFSetDirty (ObitDConCleanOTF* in, ObitImage* image) {
  in->dirty = ObitImageUnref(in->dirty);
  in->dirty = ObitImageRef(image);
}

extern ObitImage* CleanOTFGetBeam (ObitDConCleanOTF* in) {
  return ObitImageRef(in->beam);
}

extern void CleanOTFSetBeam (ObitDConCleanOTF* in, ObitImage* image) {
  in->beam = ObitImageUnref(in->beam);
  in->beam = ObitImageRef(image);
}

extern ObitImage* CleanOTFGetClean (ObitDConCleanOTF* in) {
  return ObitImageRef(in->clean);
}

extern void CleanOTFSetClean (ObitDConCleanOTF* in, ObitImage* image) {
  in->clean = ObitImageUnref(in->clean);
  in->clean = ObitImageRef(image);
}

extern void CleanOTFDefWindow (ObitDConCleanOTF* in, ObitErr *err) {
 ObitDConCleanDefWindow((ObitDConClean*)in, err);
}

extern ObitDConCleanWindow* CleanOTFGetWindow (ObitDConCleanOTF* in) {
  return ObitDConCleanWindowRef(in->window);
}

extern void CleanOTFSetWindow (ObitDConCleanOTF* in, ObitDConCleanWindow *window) {
  in->mosaic = ObitDConCleanWindowUnref(in->window);  /* Out with the old */
  in->mosaic = ObitDConCleanWindowRef(window);        /* In with the new */
}

// if (win[0]<0 then this is a round window then win[1]=radius, [2,3] = center
// else rectangular and blc=(win[0],win[1]), trc= blc=(win[2],win[3])
extern void CleanOTFAddWindow (ObitDConCleanOTF* in, int *win, 
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
}  // end CleanOTFAddWindow 

extern ObitInfoList* CleanOTFGetList (ObitDConCleanOTF* in) {
  return ObitInfoListRef(in->info);
}

extern long CleanOTFGetNiter (ObitDConCleanOTF* in) {
  return in->niter;
}

extern float CleanOTFGetGain (ObitDConCleanOTF* in) {
  return in->gain[0];
}

extern float CleanOTFGetFlux (ObitDConCleanOTF* in) {
  return in->minFlux[0];
}

extern float CleanOTFGetCleanSize (ObitDConCleanOTF* in) {
  return in->cleanSize;
}

extern int CleanOTFIsA (ObitDConCleanOTF* in) {
  return ObitDConCleanOTFIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitDConCleanOTF *me;
} CleanOTF;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitDConCleanOTF *me;
} CleanOTF;

%addmethods CleanOTF { 
  CleanOTF(char* name, 
	   ObitImage *dirty, ObitImage *beam,
	   ObitImage *clean, ObitErr *err) {
     CleanOTF *out;
     // May not be passed a beam, name of beam will be "NoBeam"
     ObitImage *lbeam=NULL;
     if (beam && strncmp(beam->name,"NoBeam",6)) lbeam = beam;
     out = (CleanOTF *) malloc(sizeof(CleanOTF));
     if(name) out->me = CleanOTFCreate(name, dirty, lbeam, clean, err);
     else  out->me = NULL;
     return out;
   }
  ~CleanOTF() {
    self->me = CleanOTFUnref(self->me);
    free(self);
  }
};



