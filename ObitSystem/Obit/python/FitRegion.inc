/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for FitRegion type                         */
/*                                                                    */
/*;  Copyright (C) 2007-2019                                          */
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
#include "ObitFitRegion.h"
#include "ObitFitModel.h"
%}

/* Definitions for Python Shadow class */
/* A copy of the struct for c */
%{
typedef struct {
  ObitFitRegion *me;
} FitRegion;
%}
/* and a copy of the struct for swig */
typedef struct {
  ObitFitRegion *me;
} FitRegion;

%inline %{
extern ObitFitRegion* newFitRegion (char* name) {
  return newObitFitRegion (name);
} // end  newFitRegion

extern ObitFitRegion* FitRegionUnref (ObitFitRegion* in) {
  if (!ObitFitRegionIsA(in)) return NULL;
  return ObitFitRegionUnref(in);
}

extern ObitFitRegion*  FitRegionRef (ObitFitRegion* in) {
  return ObitFitRegionRef(in);
}

extern ObitFitRegion* FitRegionCreate (char *name, olong corner[2], olong dim[2],
		     ofloat peak, ofloat peakResid, ofloat RMSResid,
		     ofloat fluxResid) {
  return ObitFitRegionCreate((gchar*)name, corner, dim,
		             peak, peakResid, RMSResid,
		             fluxResid, 0, NULL);
}

void CreateFitRegion(PyObject *self, char* name, long *corner, long *dim, float peak, float peakResid, 
                     float RMSResid, float fluxResid) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_FitRegion, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"CreateFitRegion: could not recover c struct");
	return;
   }
   olong lcorner[2], ldim[2];
   lcorner[0] = (olong)corner[0]; lcorner[1] = (olong)corner[1]; 
   ldim[0] = (olong)dim[0]; ldim[1] = (olong)dim[1]; 
   if (strcmp(name, "None")) ((FitRegion*)ptr)->me = FitRegionCreate(name, lcorner, ldim, (ofloat)peak, 
                                 (ofloat)peakResid, (ofloat)RMSResid, (ofloat)fluxResid);
   else  ((FitRegion*)ptr)->me = NULL;
   return;
 } // end CreateFitRegion
void DeleteFitRegion(PyObject *self) {
   void *ptr;
   int ret, flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_FitRegion, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"DeleteFitRegion: could not recover c struct");
	return;
   }
   ((FitRegion*)ptr)->me = FitRegionUnref(((FitRegion*)ptr)->me);
}// end DeleteFitRegion

ObitFitRegion* FitRegion_Get_me(PyObject *self) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_FitRegion, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"FitRegion_Get_me: could not recover c struct");
	return NULL;
   }
   return ((FitRegion*)ptr)->me;
}  // end FitRegion_Get_me

void FitRegion_Set_me(PyObject *self, ObitFitRegion *me) {
   void *ptr;
   int ret,flags=0;
   ret = SWIG_Python_ConvertPtr(self, &ptr, SWIGTYPE_p_FitRegion, flags);
   if (!SWIG_IsOK(ret)) {
	PyErr_SetString(PyExc_RuntimeError,"FitRegion_Set_me: could not recover c struct");
	return;
   }
   ((FitRegion*)ptr)->me =  ObitFitRegionRef(me);
}  // end FitRegion_Set_me

extern ObitFitRegion* FitRegionCopy  (ObitFitRegion *in, ObitFitRegion *out, 
				    ObitErr *err) {
  return ObitFitRegionCopy (in, out, err);
} // end  FitRegionCopy

extern void FitRegionSetCorner(ObitFitRegion* in, long *value){
  in->corner[0] = (olong)value[0];
  in->corner[1] = (olong)value[1];
}  // end FitRegionSetCorner

extern void FitRegionSetDim(ObitFitRegion* in, long *value){
  in->dim[0] = (olong)value[0];
  in->dim[1] = (olong)value[1];
}  // end FitRegionSetDim

extern void FitRegionSetPeak(ObitFitRegion* in, float value){
  in->peak = (ofloat)value;
}  // end FitRegionSetPeak

extern void FitRegionSetPeakResid(ObitFitRegion* in, float value){
  in->peakResid = (ofloat)value;
}  // end FitRegionSetPeakResid

extern void FitRegionSetRMSResid(ObitFitRegion* in, float value){
  in->RMSResid = (ofloat)value;
}  // end FitRegionSetRMSResid

extern void FitRegionSetFluxResid(ObitFitRegion* in, float value){
  in->fluxResid = (ofloat)value;
}  // end FitRegionSetFluxResid

extern void FitRegionSetNmodel(ObitFitRegion* in, long value){
  olong i;
  in->nmodel = (olong)value;
  // resize models
  if (in->models) g_free(in->models);
  in->models = g_malloc0(in->nmodel*sizeof(ObitFitModel*));
  for (i=0; i<in->nmodel; i++) in->models[i] = NULL;
}  // end FitRegionSetNmodel

extern void FitRegionSetModels(ObitFitRegion* in, ObitFitModel *value, long i){
  if (i>=in->nmodel) {
 	PyErr_SetString(PyExc_RuntimeError,"Attempt to set invalid model no.");
        return;
  }
  in->models[i] = ObitFitModelRef(value);
}  // end FitRegionSetModels

extern PyObject*  FitRegionGetCorner(ObitFitRegion* in){
  PyObject *outList, *o;
  int i;

  outList = PyList_New(2);
  for (i=0; i<2; i++) {
    o = PyLong_FromLong((long)in->corner[i]);
    PyList_SetItem(outList, i, o);
  }
  return outList;
}  // end FitRegionGetCorner

extern PyObject*  FitRegionGetDim(ObitFitRegion* in){
  PyObject *outList, *o;
  int i;

  outList = PyList_New(2);
  for (i=0; i<2; i++) {
    o = PyLong_FromLong((long)in->dim[i]);
    PyList_SetItem(outList, i, o);
  }
  return outList;
}  // end FitRegionGetDim

extern float FitRegionGetPeak(ObitFitRegion* in){
  return (float)in->peak;
}  // end FitRegionGetPeak

extern float FitRegionGetPeakResid(ObitFitRegion* in){
  return (float)in->peakResid;
}  // end FitRegionGetPeakResid

extern float FitRegionGetRMSResid(ObitFitRegion* in){
  return (float)in->RMSResid;
}  // end FitRegionGetRMSResid

extern float FitRegionGetFluxResid(ObitFitRegion* in){
  return (float)in->fluxResid;
}  // end FitRegionGetFluxResid

extern long FitRegionGetNmodel(ObitFitRegion* in){
  return (long)in->nmodel;
}  // end FitRegionGetNmodel

extern ObitFitModel* FitRegionGetModels(ObitFitRegion* in, int i) {
  return ObitFitModelRef(in->models[i]);
}  // end FitRegionGetModels

extern char* FitRegionGetName (ObitFitRegion* in) {
  if (ObitFitRegionIsA(in)) {
    return in->name;
  } else {
    return NULL;
  }
}

extern long FitRegionIsA (ObitFitRegion* in) {
  return (long)ObitFitRegionIsA(in);
}

%}

