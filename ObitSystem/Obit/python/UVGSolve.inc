/* $Id: UVGSolve.inc,v 1.2 2006/03/13 16:55:30 bcotton Exp $ */  
/*--------------------------------------------------------------------*/
/* Swig module description for UV data self calibration utilities     */
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
#include "ObitUVGSolve.h"
#include "ObitUVSoln.h"
%}


%inline %{
extern ObitUVGSolve* newUVGSolve (char* name) {
  return newObitUVGSolve (name);
} // end  newUVGSolve

extern ObitUVGSolve* UVGSolveCreate (char *name) {
 return ObitUVGSolveCreate(name);
}

extern ObitUVGSolve* UVGSolveCopy  (ObitUVGSolve *in, ObitUVGSolve *out, 
				    ObitErr *err) {
  return ObitUVGSolveCopy (in, out, err);
} // end  UVGSolveCopy

extern ObitUVGSolve* UVGSolveUnref (ObitUVGSolve* in) {
  if (!ObitUVGSolveIsA(in)) return NULL;
  return ObitUVGSolveUnref(in);
}

extern ObitUVGSolve*  UVGSolveRef (ObitUVGSolve* in) {
  return ObitUVGSolveRef(in);
}

extern ObitInfoList* UVGSolveGetList (ObitUVGSolve* in) {
  return ObitInfoListRef(in->info);
}

extern int UVGSolveIsA (ObitUVGSolve* in) {
  return ObitUVGSolveIsA(in);
}

extern ObitTable* UVGSolveCal (ObitUVGSolve *in, ObitUV *inUV, ObitUV *outUV, 
                               ObitUVSel *sel, ObitErr *err)
{
  return (ObitTable*)ObitUVGSolveCal (in, inUV, outUV, sel, err);
} // end UVGSolveCal


%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitUVGSolve *me;
} UVGSolve;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitUVGSolve *me;
} UVGSolve;

%addmethods UVGSolve { 
  UVGSolve(char* name) {
     UVGSolve *out;
     out = (UVGSolve *) malloc(sizeof(UVGSolve));
     if (strcmp(name, "None")) out->me = newUVGSolve(name);
     else out->me = NULL;
     return out;
   }
  ~UVGSolve() {
   if (self->me->ReferenceCount>0) 
      self->me = UVGSolveUnref(self->me);
   free(self);
  }
};

