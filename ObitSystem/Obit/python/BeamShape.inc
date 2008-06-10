/* $Id: BeamShape.inc,v 1.1 2008/05/06 13:20:14 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for BeamShape type                          */
/*                                                                    */
/*;  Copyright (C) 2008                                               */
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
#include "ObitBeamShape.h"
#include "ObitImage.h"
%}


%inline %{
extern ObitBeamShape* newBeamShape (char* name) {
  return newObitBeamShape (name);
} // end  newBeamShape

extern ObitBeamShape* BeamShapeCopy  (ObitBeamShape *in, ObitBeamShape *out, 
				    ObitErr *err) {
  return ObitBeamShapeCopy (in, out, err);
} // end  BeamShapeCopy

extern ObitBeamShape* BeamShapeUnref (ObitBeamShape* in) {
  if (!ObitBeamShapeIsA(in)) return NULL;
  return ObitBeamShapeUnref(in);
}

extern ObitBeamShape*  BeamShapeRef (ObitBeamShape* in) {
  return ObitBeamShapeRef(in);
}

extern ObitBeamShape* BeamShapeCreate (char *name, ObitImage *image, 
				       float pbmin, float antSize, 
				       int doGain) {
  gboolean ldoGain = (doGain!=0);
  return ObitBeamShapeCreate((gchar*)name, image, (ofloat)pbmin, (ofloat)antSize,
			     ldoGain);
} // end BeamShapeCreate

extern float BeamShapeGain (ObitBeamShape* in, double ra, double dec,
		            float parAng ) {
  return ObitBeamShapeGain(in, (odouble)ra, (odouble)dec, (ofloat)parAng);
} // end BeamShapeGain

extern float BeamShapeGainSym (ObitBeamShape* in, double Angle) {
  return ObitBeamShapeGainSym(in, (odouble)Angle);
} // end BeamShapeGainSym

extern char* BeamShapeGetName (ObitBeamShape* in) {
  if (ObitBeamShapeIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int BeamShapeIsA (ObitBeamShape* in) {
  return ObitBeamShapeIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitBeamShape *me;
} BeamShape;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitBeamShape *me;
} BeamShape;

%addmethods BeamShape { 
  BeamShape(char* name, ObitImage *image, float pbmin, float antSize,
            int doGain) {
     BeamShape *out;
     gboolean ldoGain = (doGain!=0);
     out = (BeamShape *) malloc(sizeof(BeamShape));
     if (strcmp(name, "None")) out->me = BeamShapeCreate((gchar*)name, image, 
                                                         (ofloat)pbmin, (ofloat)antSize,
                                			  ldoGain);
     else out->me = NULL;
     return out;
   }
  ~BeamShape() {
   if (!self) return;  // Not defined
   if (self && self->me && self->me->ReferenceCount>0) {
      self->me = BeamShapeUnref(self->me);
      free(self);
   }
  }
};

