/* $Id$ */  
/*--------------------------------------------------------------------*/
/* Swig module description for ImageMosaic type                       */
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
#include "ObitUVImager.h"
%}


%inline %{
extern ObitUVImager* newUVImager (char* name) {
  return newObitUVImager (name);
} // end  newUVImager

extern ObitUVImager* UVImagerCopy  (ObitUVImager *in, ObitUVImager *out, 
				    ObitErr *err) {
  return ObitUVImagerCopy (in, out, err);
} // end  UVImagerCopy

extern ObitUVImager* UVImagerUnref (ObitUVImager* in) {
  if (!ObitUVImagerIsA(in)) return NULL;
  return ObitUVImagerUnref(in);
}

extern ObitUVImager*  UVImagerRef (ObitUVImager* in) {
  return ObitUVImagerRef(in);
}

extern ObitUV* UVImagerGetUV (ObitUVImager* in) {
  return ObitUVRef(in->uvdata);
}

extern ObitImageMosaic* UVImagerGetMosaic (ObitUVImager* in, ObitErr *err) {
  return ObitUVImagerGetMosaic(in, err);
}

extern ObitUVImager* UVImagerCreate (char *name, ObitUV *uvData, ObitErr *err) {
 return ObitUVImagerCreate(name, uvData, err);
}

extern void UVImagerWeight (ObitUVImager* in, ObitErr *err) {
 ObitUVImagerWeight(in, err);
}

extern void UVImagerImage (ObitUVImager* in, int field, int doWeight, int doBeam, 
                            int doFlatten, ObitErr *err) {
 gboolean LdoWeight=doWeight, LdoBeam=doBeam, LdoFlatten=doFlatten;
 olong Lfield = field;
 ObitUVImagerImage(in, Lfield, LdoWeight, LdoBeam, LdoFlatten, err);
}

extern void UVImagerFlatten (ObitUVImager* in, ObitErr *err) {
 ObitUVImagerFlatten(in, err);
}

extern char* UVImagerGetName (ObitUVImager* in) {
  if (ObitUVImagerIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int UVImagerIsA (ObitUVImager* in) {
  return ObitUVImagerIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitUVImager *me;
} UVImager;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitUVImager *me;
} UVImager;

%addmethods UVImager { 
  UVImager(char* name, ObitUV *uvData, ObitErr *err) {
     UVImager *out;
     out = (UVImager *) malloc(sizeof(UVImager));
     if (strcmp(name, "None")) out->me = UVImagerCreate(name, uvData, err);
     else out->me = NULL;
     return out;
   }
  ~UVImager() {
   if (self->me->ReferenceCount>0) 
      self->me = UVImagerUnref(self->me);
   free(self);
  }
};

