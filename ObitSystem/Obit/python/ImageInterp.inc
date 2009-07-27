/* $Id: $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageInterp type                       */
/*                                                                    */
/*;  Copyright (C) 2009                                               */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

%{
#include "ObitImageInterp.h"
#include "ObitImage.h"
%}


%inline %{
extern ObitImageInterp* newImageInterp (char* name) {
  return newObitImageInterp (name);
} // end  newImageInterp

extern ObitImageInterp* ImageInterpCopy  (ObitImageInterp *in, ObitImageInterp *out, 
				    ObitErr *err) {
  return ObitImageInterpCopy (in, out, err);
} // end  ImageInterpCopy

extern ObitImageInterp* ImageInterpUnref (ObitImageInterp* in) {
  if (!ObitImageInterpIsA(in)) return NULL;
  return ObitImageInterpUnref(in);
}

extern ObitImageInterp*  ImageInterpRef (ObitImageInterp* in) {
  return ObitImageInterpRef(in);
}

extern ObitImageInterp* ImageInterpCreate (char *name, ObitImage *image, 
				     ObitErr *err) {
  return ObitImageInterpCreate((gchar*)name, image, err);
} // end ImageInterpCreate

extern float ImageInterpValue (ObitImageInterp* in, double ra, double dec,
		            float parAng, int plane, ObitErr *err) {
  return ObitImageInterpValue(in, (odouble)ra, (odouble)dec, (ofloat)parAng, (olong)plane, err);
} // end ImageInterpValue

extern long ImageInterpFindPlane (ObitImageInterp* in, double freq) {
  return (long)ObitImageInterpFindPlane(in, (odouble)freq);
} // end ImageInterpFindPlane

extern char* ImageInterpGetName (ObitImageInterp* in) {
  if (ObitImageInterpIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int ImageInterpIsA (ObitImageInterp* in) {
  return ObitImageInterpIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitImageInterp *me;
} ImageInterp;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitImageInterp *me;
} ImageInterp;

%addmethods ImageInterp { 
  ImageInterp(char* name, ObitImage *image, ObitErr *err) {
     ImageInterp *out;
     out = (ImageInterp *) malloc(sizeof(ImageInterp));
     if (strcmp(name, "None")) out->me = ImageInterpCreate((gchar*)name, image, err);
     else out->me = NULL;
     return out;
   }
  ~ImageInterp() {
   if (!self) return;  // Not defined
   if (self && self->me && self->me->ReferenceCount>0) {
      self->me = ImageInterpUnref(self->me);
      free(self);
   }
  }
};

