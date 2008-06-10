/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for FitModel type                          */
/*                                                                    */
/*;  Copyright (C) 2007                                               */
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
#include "ObitFitModel.h"
#include "ObitImageFitData.h"
#include "ObitImageDesc.h"
%}


%inline %{
extern ObitFitModel* newFitModel (char* name) {
  return newObitFitModel (name);
} // end  newFitModel

extern ObitFitModel* FitModelCopy  (ObitFitModel *in, ObitFitModel *out, 
				    ObitErr *err) {
  return ObitFitModelCopy (in, out, err);
} // end  FitModelCopy

extern ObitFitModel* FitModelUnref (ObitFitModel* in) {
  if (!ObitFitModelIsA(in)) return NULL;
  return ObitFitModelUnref(in);
}

extern ObitFitModel*  FitModelRef (ObitFitModel* in) {
  return ObitFitModelRef(in);
}

extern ObitFitModel* FitModelCreate (char *name, int type, float Peak, float DeltaX, float DeltaY, 
			             int nparm, float *parms) {
  return ObitFitModelCreate((gchar*)name, (olong)type, (ofloat)Peak, (ofloat)DeltaX, (ofloat)DeltaY, 
			   (olong)nparm, (ofloat*)parms);
}

extern int FitModelGetType(ObitFitModel* in) {
  return (int)in->type;
}  // end FitModelGetType

extern float FitModelGetPeak(ObitFitModel* in) {
  return (float)in->Peak;
}  // end FitModelGetPeak(

extern float FitModelGetDeltaX(ObitFitModel* in) {
  return (float)in->DeltaX;
}  // end FitModelGetDeltaX

extern float FitModelGetDeltaY(ObitFitModel* in) {
  return (float)in->DeltaY;
}  // end FitModelGetDeltaY

extern int FitModelGetNparm(ObitFitModel* in) {
  return (int)in->nparm;
}  // end FitModelGetNparm

extern PyObject*  FitModelGetParms(ObitFitModel* in) {
  PyObject *outList, *o;
  int i;

  outList = PyList_New(in->nparm);
  for (i=0; i<in->nparm; i++) {
    o = PyFloat_FromDouble((double)in->parms[i]);
    PyList_SetItem(outList, i, o);
  }
  return outList;
}  // end FitModelGetParms

extern float FitModelGetePeak(ObitFitModel* in) {
  return (float)in->ePeak;
}  // end FitModelGetePeak

extern float FitModelGeteDeltaX(ObitFitModel* in) {
  return (float)in->eDeltaX;
}  // end FitModelGeteDeltaX

extern float FitModelGeteDeltaY(ObitFitModel* in) {
  return (float)in->eDeltaY;
}  // end FitModelGeteDeltaY

extern PyObject* FitModelGeteParms(ObitFitModel* in) {
  PyObject *outList, *o;
  int i;

  outList = PyList_New(in->nparm);
  for (i=0; i<in->nparm; i++) {
    o = PyFloat_FromDouble((double)in->eparms[i]);
    PyList_SetItem(outList, i, o);
  }
  return outList;
}  // end FitModelGeteParms

extern void FitModelSetType(ObitFitModel* in, int value) {
  in->type = (ObitFitModelCompType)value;
}  // end FitModelSetType

extern void FitModelSetPeak(ObitFitModel* in, float value) {
  in->Peak = (ofloat)value;
}  // end FitModelSetPeak

extern void FitModelSetDeltaX(ObitFitModel* in, float value) {
  in->DeltaX = (ofloat)value;
}  // end FitModelSetDeltaX

extern void FitModelSetDeltaY(ObitFitModel* in, float value) {
  in->DeltaY = (ofloat)value;
}  // end FitModelSetDeltaY

extern void FitModelSetNparm(ObitFitModel* in, int value) {
  in->nparm = (olong)value;
  if (in->parms)  g_free(in->parms);
  if (in->eparms) g_free(in->eparms);
  if (in->nparm>0) {
    in->parms  = g_malloc(in->nparm*sizeof(ofloat));
    in->eparms = g_malloc(in->nparm*sizeof(ofloat));
  } else {
    in->parms  = NULL;
    in->eparms = NULL;
  }
}  // end FitModelSetNparm

extern void FitModelSetParms(ObitFitModel* in, float *value) {
  int i;

  for (i=0; i<in->nparm; i++) {
    in->parms[i] = (ofloat)value[i];
  }
}  // end FitModelSetParms

extern void FitModelSetePeak(ObitFitModel* in, float value) {
  in->ePeak = (ofloat)value;
}  // end FitModelSetePeak

extern void FitModelSeteDeltaX(ObitFitModel* in, float value) {
  in->eDeltaX = (ofloat)value;
}  // end FitModelSeteDeltaX

extern void FitModelSeteDeltaY(ObitFitModel* in, float value) {
  in->eDeltaY = (ofloat)value;
}  // end FitModelSeteDeltaY

extern void FitModelSeteParms(ObitFitModel* in, float *value) {
  int i;

  for (i=0; i<in->nparm; i++) {
    in->eparms[i] = (ofloat)value[i];
  }
}  // end FitModelSeteParms

extern ObitFitModel* DeconGau (ObitFitModel* in, ObitFitModel* out, ObitImageDesc *desc) {
   ofloat dgau[3][3], emaj, emin, epa, xcell;
   olong i, ret;

   /* Check */
    if ((in->type!=OBIT_FitModel_GaussMod) || (in->nparm!=3)) {
 	PyErr_SetString(PyExc_TypeError,"Input model not a Gaussian");
        return NULL;
    }
   xcell = fabs (desc->cdelt[1]);	
   ret = ObitFitModelDeconGau ((ofloat)in->parms[0],   (ofloat)in->parms[1],   (ofloat)(in->parms[2]*RAD2DG), 
			       (ofloat)in->eparms[0],  (ofloat)in->eparms[1],  (ofloat)(in->eparms[2]*RAD2DG),
			       (ofloat)desc->beamMaj/xcell,  (ofloat)desc->beamMin/xcell,  
	                       (ofloat)desc->beamPA, dgau);
  if (ret>1) {
 	//PyErr_SetString(PyExc_RuntimeError,"Deconvolution failed");
	 out->type  = -10;
        return out;
  }
   dgau[0][2]  *= DG2RAD;
   dgau[1][2]  *= DG2RAD;
   dgau[2][2]  *= DG2RAD;
   emaj = fabs(dgau[1][0] - dgau[2][0]) * 0.5;
   emin = fabs(dgau[1][1] - dgau[2][1]) * 0.5;
   epa  = fabs(dgau[1][2] - dgau[2][2]) * 0.5;
   if (epa>G_PI)  epa -= 2.0*G_PI;
   if (epa<-G_PI) epa += 2.0*G_PI;
   out->type    = in->type;
   out->nparm   = in->nparm;
   out->Peak    = in->Peak;
   out->DeltaX  = in->DeltaX;
   out->DeltaY  = in->DeltaY;
   out->ePeak   = in->ePeak;
   out->eDeltaX = in->eDeltaX;
   out->eDeltaY = in->eDeltaY;
   for (i=0; i<3; i++) out->parms[i] = dgau[0][i];
   out->eparms[0] = emaj;
   out->eparms[1] = emin;
   out->eparms[2] = epa;

   return out;
} // end DeconGau

// Determine errors for Gaussian model
extern void FitModelGaussErr (ObitFitModel* in, ObitImageDesc *desc, float irms) {
   ofloat beam[3], xcell;

   /* Check that Gaussian */
    if ((in->type!=OBIT_FitModel_GaussMod) || (in->nparm!=3)) {
 	PyErr_SetString(PyExc_TypeError,"Input model not a Gaussian");
        return;
    }
   xcell = fabs (desc->cdelt[1]);	
   beam[0] = (ofloat)desc->beamMaj/xcell;
   beam[1] = (ofloat)desc->beamMin/xcell;
   beam[2] = (ofloat)desc->beamPA;

   ObitImageFitDataGaussErr ((ofloat)in->Peak,
                             (ofloat)in->parms[0], (ofloat)in->parms[1], (ofloat)(in->parms[2]), 
			     (ofloat)irms,  (ofloat*)&beam[0],
			     (ofloat*)&in->ePeak, (ofloat*)&in->eDeltaX,  (ofloat*)&in->eDeltaY, 
			     (ofloat*)&in->eparms[0], (ofloat*)&in->eparms[1],  (ofloat*)&in->eparms[2]);
   return;
} // end FitModelGaussErr

extern char* FitModelGetName (ObitFitModel* in) {
  if (ObitFitModelIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern int FitModelIsA (ObitFitModel* in) {
  return ObitFitModelIsA(in);
}

%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitFitModel *me;
} FitModel;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitFitModel *me;
} FitModel;

%addmethods FitModel { 
  FitModel(char* name, int type, float Peak, float DeltaX, float DeltaY, int nparm, float *parms) {
     FitModel *out;
     out = (FitModel *) malloc(sizeof(FitModel));
     if (strcmp(name, "None")) out->me = FitModelCreate(name, type, Peak, DeltaX, DeltaY, 
			                                nparm, parms);
     else out->me = NULL;
     return out;
   }
  ~FitModel() {
   if (!self) return;  // Not defined
   if (self && self->me && self->me->ReferenceCount>0) {
      self->me = FitModelUnref(self->me);
      free(self);
   }
  }
};

