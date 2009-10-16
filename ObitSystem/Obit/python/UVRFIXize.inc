/* $Id:  $ */  
/*--------------------------------------------------------------------*/
/* Swig module description for UV data self calibration utilities     */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

%{
#include "ObitUVRFIXize.h"
%}


%inline %{
extern ObitUVRFIXize* newUVRFIXize (char* name) {
  return newObitUVRFIXize (name);
} // end  newUVRFIXize

extern ObitUVRFIXize* UVRFIXizeCreate (char *name, ObitUV *inUV,
	ObitUV *residUV, ObitUV *outUV) {
 return ObitUVRFIXizeCreate(name, inUV, residUV, outUV);
}

extern ObitUVRFIXize* UVRFIXizeUnref (ObitUVRFIXize* in) {
  if (!ObitUVRFIXizeIsA(in)) return NULL;
  return ObitUVRFIXizeUnref(in);
}

extern ObitUVRFIXize*  UVRFIXizeRef (ObitUVRFIXize* in) {
  return ObitUVRFIXizeRef(in);
}

extern void UVRFIXizeCounterRot (ObitUVRFIXize* in, ObitErr *err) {
  ObitUVRFIXizeCounterRot(in, err);
}

extern void UVRFIXizeFilter (ObitUVRFIXize* in, ObitErr *err) {
  ObitUVRFIXizeFilter(in, err);
}

extern void UVRFIXizeCorrect (ObitUVRFIXize* in, ObitErr *err) {
  ObitUVRFIXizeCorrect(in, err);
}

extern ObitInfoList* UVRFIXizeGetList (ObitUVRFIXize* in) {
  return ObitInfoListRef(in->info);
}

extern ObitSkyModel* UVRFIXizeGetRFI (ObitUVRFIXize* in) {
  return ObitUVRef(in->RFIUV);
}

extern int UVRFIXizeIsA (ObitUVRFIXize* in) {
  return ObitUVRFIXizeIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitUVRFIXize *me;
} UVRFIXize;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitUVRFIXize *me;
} UVRFIXize;

%addmethods UVRFIXize { 
  UVRFIXize(char* name) {
     UVRFIXize *out;
     out = (UVRFIXize *) malloc(sizeof(UVRFIXize));
     if (strcmp(name, "None")) out->me = newUVRFIXize(name);
     else out->me = NULL;
     return out;
   }
  ~UVRFIXize() {
   if (self->me->ReferenceCount>0) 
      self->me = UVRFIXizeUnref(self->me);
   free(self);
  }
};

