/* $Id$ */  
/*--------------------------------------------------------------------*/
/* Swig module description for UV data self calibration utilities     */
/*                                                                    */
/*;  Copyright (C) 2005-2008                                          */
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
#include "ObitUVSelfCal.h"
#include "ObitUVSoln.h"
%}


%inline %{
extern ObitUVSelfCal* newUVSelfCal (char* name) {
  return newObitUVSelfCal (name);
} // end  newUVSelfCal

extern ObitUVSelfCal* UVSelfCalCreate (char *name, ObitSkyModel *skyModel) {
 return ObitUVSelfCalCreate(name, skyModel);
}

extern ObitUVSelfCal* UVSelfCalCopy  (ObitUVSelfCal *in, ObitUVSelfCal *out, 
				    ObitErr *err) {
  return ObitUVSelfCalCopy (in, out, err);
} // end  UVSelfCalCopy

extern ObitUVSelfCal* UVSelfCalUnref (ObitUVSelfCal* in) {
  if (!ObitUVSelfCalIsA(in)) return NULL;
  return ObitUVSelfCalUnref(in);
}

extern ObitUVSelfCal*  UVSelfCalRef (ObitUVSelfCal* in) {
  return ObitUVSelfCalRef(in);
}

extern ObitInfoList* UVSelfCalGetList (ObitUVSelfCal* in) {
  return ObitInfoListRef(in->info);
}

extern ObitSkyModel* UVSelfCalGetSkyModel (ObitUVSelfCal* in) {
  return ObitSkyModelRef(in->skyModel);
}

extern void UVSelfCalSetSkyModel (ObitUVSelfCal* in, ObitSkyModel *skyModel, 
                                     ObitErr *err) {
  in->skyModel = ObitSkyModelUnref(in->skyModel);  /* Out with the old */
  in->skyModel = ObitSkyModelRef(skyModel);        /* In with the new */
}

extern int UVSelfCalIsA (ObitUVSelfCal* in) {
  return ObitUVSelfCalIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitUVSelfCal *me;
} UVSelfCal;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitUVSelfCal *me;
} UVSelfCal;

%addmethods UVSelfCal { 
  UVSelfCal(char* name) {
     UVSelfCal *out;
     out = (UVSelfCal *) malloc(sizeof(UVSelfCal));
     if (strcmp(name, "None")) out->me = newUVSelfCal(name);
     else out->me = NULL;
     return out;
   }
  ~UVSelfCal() {
   if (self->me->ReferenceCount>0) 
      self->me = UVSelfCalUnref(self->me);
   free(self);
  }
};

