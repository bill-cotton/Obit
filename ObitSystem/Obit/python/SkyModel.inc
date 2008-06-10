/* $Id: SkyModel.inc,v 1.5 2007/11/06 13:00:31 bcotton Exp $    */  
/*--------------------------------------------------------------------*/
/* Swig module description for SkyModel type                          */
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
#include "ObitSkyModel.h"
%}


%inline %{
extern ObitSkyModel* newSkyModel (char* name) {
  return newObitSkyModel (name);
} // end  newSkyModel

extern ObitSkyModel* SkyModelCopy  (ObitSkyModel *in, ObitSkyModel *out, 
				    ObitErr *err) {
  return ObitSkyModelCopy (in, out, err);
} // end  SkyModelCopy

extern ObitSkyModel* SkyModelUnref (ObitSkyModel* in) {
  if (!ObitSkyModelIsA(in)) return NULL;
  return ObitSkyModelUnref(in);
}

extern ObitSkyModel*  SkyModelRef (ObitSkyModel* in) {
  return ObitSkyModelRef(in);
}

extern ObitInfoList* SkyModelGetList (ObitSkyModel* in) {
  return ObitInfoListRef(in->info);
}

extern ObitImageMosaic* SkyModelGetImageMosaic (ObitSkyModel* in) {
  return ObitImageMosaicRef(in->mosaic);
}

extern void SkyModelSetImageMosaic (ObitSkyModel* in, ObitImageMosaic *mosaic, 
                              ObitErr *err) {
  in->mosaic = ObitImageMosaicUnref(in->mosaic);  /* Out with the old */
  in->mosaic = ObitImageMosaicRef(mosaic);        /* In with the new */
}

extern ObitSkyModel* SkyModelCreate (char *name, ObitImageMosaic* mosaic) {
 return ObitSkyModelCreate(name, mosaic);
}

extern void SkyModelSubUV (ObitSkyModel* in, ObitUV *inData, ObitUV *outData, ObitErr *err) {
 ObitSkyModelSubUV(in, inData, outData, err);
}

extern void SkyModelDivUV (ObitSkyModel* in, ObitUV *inData, ObitUV *outData, ObitErr *err) {
 ObitSkyModelDivUV(in, inData, outData, err);
}

extern void SkyModelCompressCC (ObitSkyModel* in, ObitErr *err) {
   ObitSkyModelCompressCC (in, err);
}

extern char* SkyModelGetName (ObitSkyModel* in) {
  if (ObitSkyModelIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int SkyModelIsA (ObitSkyModel* in) {
  return ObitSkyModelIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitSkyModel *me;
} SkyModel;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitSkyModel *me;
} SkyModel;

%addmethods SkyModel { 
  SkyModel(char* name) {
     SkyModel *out;
     out = (SkyModel *) malloc(sizeof(SkyModel));
     if (strcmp(name, "None")) out->me = newSkyModel(name);
     else out->me = NULL;
     return out;
   }
  ~SkyModel() {
   if (self->me-> ReferenceCount>0) 
      self->me = SkyModelUnref(self->me);
   free(self);
  }
};

