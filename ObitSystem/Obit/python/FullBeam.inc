/* $Id: $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for FullBeam type                          */
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
#include "ObitFullBeam.h"
#include "ObitImage.h"
%}


%inline %{
extern ObitFullBeam* newFullBeam (char* name) {
  return newObitFullBeam (name);
} // end  newFullBeam

extern ObitFullBeam* FullBeamCopy  (ObitFullBeam *in, ObitFullBeam *out, 
				    ObitErr *err) {
  return ObitFullBeamCopy (in, out, err);
} // end  FullBeamCopy

extern ObitFullBeam* FullBeamUnref (ObitFullBeam* in) {
  if (!ObitFullBeamIsA(in)) return NULL;
  return ObitFullBeamUnref(in);
}

extern ObitFullBeam*  FullBeamRef (ObitFullBeam* in) {
  return ObitFullBeamRef(in);
}

extern ObitFullBeam* FullBeamCreate (char *name, ObitImage *image, 
				     ObitErr *err) {
  return ObitFullBeamCreate((gchar*)name, NULL, image, err);
} // end FullBeamCreate

extern float FullBeamValue (ObitFullBeam* in, double dra, double ddec,
		            float parAng, int plane, ObitErr *err) {
  return ObitFullBeamValue(in, (odouble)dra, (odouble)ddec, (ofloat)parAng, (olong)plane, err);
} // end FullBeamGain

extern long FullBeamFindPlane (ObitFullBeam* in, double freq) {
  return (long)ObitFullBeamFindPlane(in, (odouble)freq);
} // end FullBeamFindPlane

extern char* FullBeamGetName (ObitFullBeam* in) {
  if (ObitFullBeamIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int FullBeamIsA (ObitFullBeam* in) {
  return ObitFullBeamIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitFullBeam *me;
} FullBeam;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitFullBeam *me;
} FullBeam;

%addmethods FullBeam { 
  FullBeam(char* name, ObitImage *image, ObitErr *err) {
     FullBeam *out;
     out = (FullBeam *) malloc(sizeof(FullBeam));
     if (strcmp(name, "None")) out->me = FullBeamCreate((gchar*)name, image, err);
     else out->me = NULL;
     return out;
   }
  ~FullBeam() {
   if (!self) return;  // Not defined
   if (self && self->me && self->me->ReferenceCount>0) {
      self->me = FullBeamUnref(self->me);
      free(self);
   }
  }
};

